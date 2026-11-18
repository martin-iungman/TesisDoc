# EXPLORATORIO - embeddings de ADN: fine-tuning de DNABERT-2 contra una CNN ajustada
#
# Corre en Colab con GPU (finetuning_colab.ipynb). Localmente solo se usa con
# --prueba (subconjunto chico, CPU) para verificar que el codigo anda.
#
# Pregunta: si DNABERT-2 ajustado de punta a punta supera a ridge sobre sus
# embeddings congelados, y si esa ganancia es del preentrenamiento o de
# cualquier modelo profundo supervisado bien entrenado (la CNN).
#
# Diseno (igual para los dos modelos, para que la comparacion sea pareja):
#   - Mismo split por cromosoma que todo lo demas (columna fold del dataset).
#     Para cada fold de test, la validacion interna (parada temprana y
#     eleccion de hiperparametros) es ~10% de los promotores de entrenamiento,
#     agrupados por cromosoma; nunca se mira el fold de test.
#   - Multitarea con perdida enmascarada, tres salidas:
#       pct_rep1, pct_rep2  percentil dentro de cada replica (todos los
#                           promotores medidos en esa replica), como
#                           mejoras.py por defecto;
#       pct_prom            promedio de los percentiles de las dos replicas,
#                           calculados en los promotores medidos en ambas,
#                           como mejoras.py --objetivo promedio.
#     Diferencia con ridge: ahi cada objetivo se ajusta por separado (y el
#     del promedio solo con los 8889 medidos en ambas).
#   - Busqueda de hiperparametros en la validacion interna de cada fold:
#     CNN, N_CONFIGS configuraciones al azar; DNABERT-2, LRS_DNABERT2.
#     Con la mejor se entrenan SEMILLAS semillas y se promedian (_s0 = solo
#     la primera semilla, para comparar contra un unico modelo).
#   - CNN: bloques conv-BN-ReLU-maxpool-dropout y capa densa sobre el mapa
#     aplanado (conserva posicion); aumento por desplazamiento de +-d pb.
#     Sin aumento por complemento reverso: la secuencia esta orientada.
#   - DNABERT-2: modelo completo ajustado, cabezal lineal sobre el pooling
#     por las mismas tres ventanas que los embeddings congelados.
#
# Cada (modelo, fold, semilla) guarda sus predicciones del fold de test en
# <carpeta>/ft/<modelo>/; si el archivo existe no se recalcula (Colab se
# desconecta). ensamblar() escribe, en el formato de mejoras.py:
#   <carpeta>/ft/predicciones_replica/<modelo>.tsv   (seq_id, Rep 1, Rep 2)
#   <carpeta>/ft/predicciones_promedio/<modelo>.tsv  (seq_id, Promedio)
# que se copian a figures/exploratorio_embeddings_adn/mejoras/predicciones/ y
# mejoras_promedio/predicciones/ respectivamente.
#
# Prueba local (desde la raiz de TesisDoc; requiere torch y, para DNABERT-2,
# transformers==4.46.3 y einops):
#   python R/71_exploratorio_embeddings_adn/finetuning.py --prueba --carpeta <tmp> --modelos cnn

import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluar as ev  # noqa: E402

OBJETIVOS = ["pct_rep1", "pct_rep2", "pct_prom"]
ESCALA = 100 / np.sqrt(12)  # desvio de una uniforme 0-100
SEMILLAS = (0, 1, 2)
FRAC_VAL = 0.1
N_CONFIGS = 12
LRS_DNABERT2 = (2e-5, 5e-5)
VENT = list(ev.VENTANAS.values())


# --- Datos -------------------------------------------------------------------------

def objetivos(df):
    """Percentiles como en mejoras.py (dataset_percentil y dataset_promedio)."""
    Y = pd.DataFrame(index=df.index)
    Y["pct_rep1"] = 100 * df.mean_rep1.rank(pct=True)
    Y["pct_rep2"] = 100 * df.mean_rep2.rank(pct=True)
    ambas = df.mean_rep1.notna() & df.mean_rep2.notna()
    p1 = 100 * df.mean_rep1[ambas].rank(pct=True)
    p2 = 100 * df.mean_rep2[ambas].rank(pct=True)
    Y["pct_prom"] = np.nan
    Y.loc[ambas, "pct_prom"] = (p1 + p2) / 2
    return Y[OBJETIVOS]


def particion(df, f):
    """test = fold f; validacion = cromosomas al azar (fijo por fold) del resto
    hasta juntar FRAC_VAL de los promotores de entrenamiento."""
    te = (df.fold == f).to_numpy()
    resto = df.loc[~te]
    rng = np.random.default_rng(100 + int(f))
    croms = rng.permutation(sorted(resto.chromosome.unique()))
    val, n = [], 0
    for c in croms:
        if n >= FRAC_VAL * len(resto):
            break
        val.append(c)
        n += (resto.chromosome == c).sum()
    va = (~te) & df.chromosome.isin(val).to_numpy()
    tr = (~te) & ~va
    return np.where(tr)[0], np.where(va)[0], np.where(te)[0]


def fijar_semilla(s):
    random.seed(s)
    np.random.seed(s)
    torch.manual_seed(s)
    torch.cuda.manual_seed_all(s)


# --- Entrenamiento generico ----------------------------------------------------------

def perdida(pred, y):
    m = ~torch.isnan(y)
    return ((pred[m] - y[m]) ** 2).mean()


def predecir(red, lote_fn, idx, dev, tam=256):
    red.eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(idx), tam):
            with torch.autocast(device_type=dev, dtype=torch.float16, enabled=dev == "cuda"):
                out.append(red(**lote_fn(idx[i:i + tam], False)).float().cpu())
    return torch.cat(out)


def ajustar(red, lote_fn, Y, tr, va, opt, dev, epocas_max, paciencia, tam_lote, sched=None):
    """Entrena con parada temprana sobre la perdida de validacion. Y: tensor
    n x 3 escalado, con NaN donde falta el objetivo."""
    amp = dev == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    mejor, estado, sin_mejora, ep = np.inf, None, 0, 0
    for ep in range(1, epocas_max + 1):
        red.train()
        perm = np.random.permutation(tr)
        for i in range(0, len(perm), tam_lote):
            b = perm[i:i + tam_lote]
            with torch.autocast(device_type=dev, dtype=torch.float16, enabled=amp):
                out = red(**lote_fn(b, True))
            l = perdida(out.float(), Y[b].to(dev))
            opt.zero_grad(set_to_none=True)
            scaler.scale(l).backward()
            scaler.unscale_(opt)
            nn.utils.clip_grad_norm_(red.parameters(), 1.0)
            scaler.step(opt)
            scaler.update()
            if sched is not None:
                sched.step()
        lv = perdida(predecir(red, lote_fn, va, dev), Y[va]).item()
        if lv < mejor - 1e-4:
            mejor, sin_mejora = lv, 0
            estado = {k: v.detach().clone() for k, v in red.state_dict().items()}
        else:
            sin_mejora += 1
            if sin_mejora >= paciencia:
                break
    red.load_state_dict(estado)
    return mejor, ep


# --- CNN -----------------------------------------------------------------------------

class CNN(nn.Module):
    def __init__(self, filtros, k1, n_bloques, drop, denso, L=252):
        super().__init__()
        capas, c_in, L_out = [], 4, L
        for i in range(n_bloques):
            k = k1 if i == 0 else 5
            capas += [nn.Conv1d(c_in, filtros, k, padding=k // 2), nn.BatchNorm1d(filtros), nn.ReLU(),
                      nn.MaxPool1d(2), nn.Dropout(drop)]
            c_in, L_out = filtros, L_out // 2
        self.conv = nn.Sequential(*capas)
        self.cabeza = nn.Sequential(nn.Flatten(), nn.Linear(filtros * L_out, denso), nn.BatchNorm1d(denso),
                                    nn.ReLU(), nn.Dropout(drop), nn.Linear(denso, 3))

    def forward(self, x):
        return self.cabeza(self.conv(x))


ESPACIO_CNN = {"filtros": [128, 256], "k1": [7, 11, 15, 19], "n_bloques": [3, 4], "drop": [0.2, 0.3, 0.4],
               "denso": [128, 256], "lr": [1e-3, 5e-4, 2e-4], "wd": [1e-4, 1e-3, 1e-2], "desp": [0, 3]}


def configs_cnn(f, n):
    rng = np.random.default_rng(200 + int(f))
    return [{k: v[rng.integers(len(v))] for k, v in ESPACIO_CNN.items()} for _ in range(n)]


def lote_cnn(onehot, dev, desp):
    X = torch.tensor(onehot)  # n x 4 x 252

    def fn(idx, entrenando):
        x = X[idx]
        if entrenando and desp > 0:
            xp = nn.functional.pad(x, (desp, desp))
            s = np.random.randint(0, 2 * desp + 1, size=len(idx))
            x = torch.stack([xp[j, :, s[j]:s[j] + x.shape[-1]] for j in range(len(idx))])
        return {"x": x.to(dev, non_blocking=True)}
    return fn


def entrenar_cnn(cfg, onehot, Y, tr, va, dev, semilla, epocas_max, paciencia):
    fijar_semilla(semilla)
    red = CNN(cfg["filtros"], cfg["k1"], cfg["n_bloques"], cfg["drop"], cfg["denso"]).to(dev)
    opt = torch.optim.AdamW(red.parameters(), lr=cfg["lr"], weight_decay=cfg["wd"])
    fn = lote_cnn(onehot, dev, cfg["desp"])
    lv, ep = ajustar(red, fn, Y, tr, va, opt, dev, epocas_max, paciencia, tam_lote=128)
    return red, fn, lv, ep


# --- DNABERT-2 -------------------------------------------------------------------------

def spans_nucleotidicos(tok, ids):
    """Igual que en extraer_embeddings_colab.ipynb."""
    toks = tok.convert_ids_to_tokens(ids)
    keep, spans, pos = [], [], 0
    for i, t in enumerate(toks):
        t = t.replace("##", "")
        if t and set(t) <= set("ACGT"):
            keep.append(i)
            spans.append((pos, pos + len(t)))
            pos += len(t)
    assert pos == 252, f"los tokens cubren {pos} pb"
    return np.array(keep), np.array(spans)


def tokenizar(tok, seqs):
    """ids por secuencia y una matriz tokens x 3 de pesos de pooling por ventana
    (1/n para los tokens de secuencia cuyo centro cae en la ventana)."""
    ids, pesos = [], []
    for s in seqs:
        i = tok(s)["input_ids"]
        keep, spans = spans_nucleotidicos(tok, i)
        centro = spans.mean(1)
        w = np.zeros((len(i), len(VENT)), np.float32)
        for j, (a, b) in enumerate(VENT):
            sel = keep[(centro >= a) & (centro < b)]
            if len(sel):
                w[sel, j] = 1 / len(sel)
        ids.append(np.array(i))
        pesos.append(w)
    return ids, pesos


def lote_dnabert2(ids, pesos, pad_id, dev):
    def fn(idx, entrenando):
        L = max(len(ids[i]) for i in idx)
        I = torch.full((len(idx), L), pad_id, dtype=torch.long)
        M = torch.zeros((len(idx), L), dtype=torch.long)
        W = torch.zeros((len(idx), L, len(VENT)))
        for j, i in enumerate(idx):
            n = len(ids[i])
            I[j, :n] = torch.from_numpy(ids[i])
            M[j, :n] = 1
            W[j, :n] = torch.from_numpy(pesos[i])
        return {"input_ids": I.to(dev), "attention_mask": M.to(dev), "pesos": W.to(dev)}
    return fn


class DNABERT2Reg(nn.Module):
    def __init__(self, base, d=768, drop=0.1):
        super().__init__()
        self.base = base
        self.drop = nn.Dropout(drop)
        self.cabeza = nn.Linear(len(VENT) * d, 3)

    def forward(self, input_ids, attention_mask, pesos):
        h = self.base(input_ids, attention_mask=attention_mask)[0]  # B x L x d (repaddeado)
        v = torch.einsum("bld,blw->bwd", h.float(), pesos).flatten(1)
        return self.cabeza(self.drop(v))


NOMBRE_DNABERT2 = "zhihan1996/DNABERT-2-117M"


def cargar_dnabert2():
    from transformers import AutoModel, AutoTokenizer
    from transformers.models.bert.configuration_bert import BertConfig
    tok = AutoTokenizer.from_pretrained(NOMBRE_DNABERT2, trust_remote_code=True)
    cfg = BertConfig.from_pretrained(NOMBRE_DNABERT2)
    return tok, lambda: AutoModel.from_pretrained(NOMBRE_DNABERT2, trust_remote_code=True, config=cfg)


def entrenar_dnabert2(nuevo_base, fn, Y, tr, va, dev, lr, semilla, epocas_max, paciencia, tam_lote=32):
    from transformers import get_linear_schedule_with_warmup
    fijar_semilla(semilla)
    red = DNABERT2Reg(nuevo_base()).to(dev)
    grupos = [{"params": red.base.parameters(), "lr": lr},
              {"params": red.cabeza.parameters(), "lr": 1e-3}]
    opt = torch.optim.AdamW(grupos, weight_decay=0.01)
    pasos = epocas_max * int(np.ceil(len(tr) / tam_lote))
    sched = get_linear_schedule_with_warmup(opt, int(0.1 * pasos), pasos)
    lv, ep = ajustar(red, fn, Y, tr, va, opt, dev, epocas_max, paciencia, tam_lote, sched)
    return red, lv, ep


# --- Orquestacion ------------------------------------------------------------------------

def _guardar(path, pred, info=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, pred)
    if info is not None:
        path.with_suffix(".json").write_text(json.dumps(info, indent=1, default=float))


def correr_cnn(df, carpeta, dev=None, n_configs=N_CONFIGS, semillas=SEMILLAS, epocas_max=60, paciencia=8):
    dev = dev or ("cuda" if torch.cuda.is_available() else "cpu")
    base = Path(carpeta) / "ft" / "cnn_ajustada"
    Y = torch.tensor(((objetivos(df) - 50) / ESCALA).to_numpy(np.float32))
    onehot = np.eye(4, dtype=np.float32)[ev.codificar(df.seq)].transpose(0, 2, 1).copy()
    for f in sorted(df.fold.unique()):
        tr, va, te = particion(df, f)
        busq = base / f"busqueda_fold{f}.json"
        if busq.exists():
            cfg = json.loads(busq.read_text())["elegida"]
        else:
            t0, res = time.time(), []
            for c in configs_cnn(f, n_configs):
                red, fn, lv, ep = entrenar_cnn(c, onehot, Y, tr, va, dev, semillas[0], epocas_max, paciencia)
                res.append({**c, "val_mse": lv, "epocas": ep})
                if lv <= min(r["val_mse"] for r in res):  # la mejor hasta ahora: su prediccion es la semilla 0
                    _guardar(base / f"fold{f}_s{semillas[0]}.npy", predecir(red, fn, te, dev).numpy(),
                             {**c, "val_mse": lv, "epocas": ep})
                print(f"  CNN fold {f}: {c} -> val {lv:.4f} ({ep} ep)", flush=True)
            cfg = min(res, key=lambda r: r["val_mse"])
            busq.write_text(json.dumps({"elegida": cfg, "todas": res}, indent=1, default=float))
            print(f"CNN fold {f}: busqueda en {time.time() - t0:.0f} s, elegida {cfg}", flush=True)
        cfg = {k: cfg[k] for k in ESPACIO_CNN}
        for s in semillas:
            p = base / f"fold{f}_s{s}.npy"
            if p.exists():
                continue
            red, fn, lv, ep = entrenar_cnn(cfg, onehot, Y, tr, va, dev, s, epocas_max, paciencia)
            _guardar(p, predecir(red, fn, te, dev).numpy(), {**cfg, "val_mse": lv, "epocas": ep})
            print(f"CNN fold {f} semilla {s}: val {lv:.4f} ({ep} ep)", flush=True)


def correr_dnabert2(df, carpeta, dev=None, lrs=LRS_DNABERT2, semillas=SEMILLAS, epocas_max=8, paciencia=2,
                    tam_lote=32):
    dev = dev or ("cuda" if torch.cuda.is_available() else "cpu")
    base = Path(carpeta) / "ft" / "dnabert2_ft"
    Y = torch.tensor(((objetivos(df) - 50) / ESCALA).to_numpy(np.float32))
    tok, nuevo_base = cargar_dnabert2()
    ids, pesos = tokenizar(tok, df.seq.tolist())
    fn = lote_dnabert2(ids, pesos, tok.pad_token_id, dev)
    for f in sorted(df.fold.unique()):
        tr, va, te = particion(df, f)
        busq = base / f"busqueda_fold{f}.json"
        if busq.exists():
            lr = json.loads(busq.read_text())["elegida"]["lr"]
        else:
            res = []
            for lr in lrs:
                t0 = time.time()
                red, lv, ep = entrenar_dnabert2(nuevo_base, fn, Y, tr, va, dev, lr, semillas[0], epocas_max,
                                                paciencia, tam_lote)
                res.append({"lr": lr, "val_mse": lv, "epocas": ep})
                if lv <= min(r["val_mse"] for r in res):
                    _guardar(base / f"fold{f}_s{semillas[0]}.npy", predecir(red, fn, te, dev).numpy(), res[-1])
                print(f"  DNABERT-2 fold {f}: lr {lr} -> val {lv:.4f} ({ep} ep, {time.time() - t0:.0f} s)",
                      flush=True)
                del red
                torch.cuda.empty_cache()
            lr = min(res, key=lambda r: r["val_mse"])["lr"]
            busq.write_text(json.dumps({"elegida": {"lr": lr}, "todas": res}, indent=1, default=float))
        for s in semillas:
            p = base / f"fold{f}_s{s}.npy"
            if p.exists():
                continue
            t0 = time.time()
            red, lv, ep = entrenar_dnabert2(nuevo_base, fn, Y, tr, va, dev, lr, s, epocas_max, paciencia, tam_lote)
            _guardar(p, predecir(red, fn, te, dev).numpy(), {"lr": lr, "val_mse": lv, "epocas": ep})
            print(f"DNABERT-2 fold {f} semilla {s}: val {lv:.4f} ({ep} ep, {time.time() - t0:.0f} s)", flush=True)
            del red
            torch.cuda.empty_cache()


def r2(y, p):
    ok = ~np.isnan(y)
    y, p = y[ok], p[ok]
    return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)


def ensamblar(df, carpeta, modelos=("cnn_ajustada", "dnabert2_ft"), semillas=SEMILLAS):
    """Junta las predicciones por fold en tablas con el formato de mejoras.py
    (<modelo> = promedio de semillas, <modelo>_s0 = primera semilla) y devuelve
    el R2 de cada una contra cada objetivo."""
    Yo = objetivos(df)
    ambas = Yo.pct_prom.notna().to_numpy()
    raiz = Path(carpeta) / "ft"
    (raiz / "predicciones_replica").mkdir(parents=True, exist_ok=True)
    (raiz / "predicciones_promedio").mkdir(parents=True, exist_ok=True)
    filas = []
    for m in modelos:
        for sufijo, sem in ((f"_s{semillas[0]}", semillas[:1]), ("", semillas)):
            P = np.full((len(df), 3), np.nan)
            completo = True
            for f in sorted(df.fold.unique()):
                te = (df.fold == f).to_numpy()
                arch = [raiz / m / f"fold{f}_s{s}.npy" for s in sem]
                if not all(a.exists() for a in arch):
                    completo = False
                    break
                P[te] = np.mean([np.load(a) for a in arch], axis=0)
            if not completo:
                print(f"{m}{sufijo}: faltan folds o semillas, no se ensambla")
                continue
            P = P * ESCALA + 50
            nombre = m + sufijo
            pd.DataFrame({"seq_id": df.seq_id, "Rep 1": P[:, 0], "Rep 2": P[:, 1]}).to_csv(
                raiz / "predicciones_replica" / f"{nombre}.tsv", sep="\t", index=False)
            pd.DataFrame({"seq_id": df.seq_id[ambas], "Promedio": P[ambas, 2]}).to_csv(
                raiz / "predicciones_promedio" / f"{nombre}.tsv", sep="\t", index=False)
            filas.append({"modelo": nombre, **{o: r2(Yo[o].to_numpy(), P[:, j]) for j, o in enumerate(OBJETIVOS)}})
    return pd.DataFrame(filas)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--carpeta", required=True, help="donde escribir ft/")
    ap.add_argument("--modelos", default="cnn,dnabert2")
    ap.add_argument("--prueba", action="store_true", help="subconjunto chico, 1 semilla, pocas epocas")
    a = ap.parse_args()
    df = ev.cargar()
    kw_cnn, kw_bert = {}, {}
    if a.prueba:
        df = df.groupby("fold", group_keys=False).apply(lambda g: g.sample(80, random_state=0)).reset_index(drop=True)
        kw_cnn = {"n_configs": 2, "semillas": (0,), "epocas_max": 2, "paciencia": 1}
        kw_bert = {"lrs": (5e-5,), "semillas": (0,), "epocas_max": 1, "paciencia": 1, "tam_lote": 8}
    mods = a.modelos.split(",")
    if "cnn" in mods:
        correr_cnn(df, a.carpeta, **kw_cnn)
    if "dnabert2" in mods:
        correr_dnabert2(df, a.carpeta, **kw_bert)
    print(ensamblar(df, a.carpeta, semillas=kw_cnn.get("semillas", SEMILLAS)).round(3))


if __name__ == "__main__":
    main()
