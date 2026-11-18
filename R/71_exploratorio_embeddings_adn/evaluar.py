# EXPLORATORIO - embeddings de modelos de lenguaje de ADN (evaluacion)
#
# Pregunta: cuanto de la actividad media (explicable) de cada promotor puede
# predecirse desde la secuencia de 252 pb sola, y cuanto agrega una
# representacion preentrenada sobre el genoma respecto de bases simples.
#
# Todo usa el mismo split: los 5 folds por cromosoma del dataset
# (preparar_dataset.py). Cada promotor recibe una prediccion "fuera de fold"
# hecha por un modelo que no vio su cromosoma. Cada replica por separado:
# se entrena con la replica r y se evalua contra la misma replica r
# (r, Spearman, R2 fuera de fold) y, en los promotores medidos en ambas,
# contra la otra replica.
#
# Techo: con un predictor perfecto del valor "verdadero" y confiabilidad
# parecida en las dos replicas, el R2 esperado contra una replica es ~ r
# entre replicas (no r^2), y la correlacion esperada es ~ sqrt(r). Se
# reporta frac_techo = R2 / r_entre_replicas.
#
# Modelos de base:
#   aditivo            OLS con las 8 caracteristicas binarias de R/38_aditividad
#   aditivo_comp       aditivo + G+C y CpG o/e continuos
#   kmeros_1a5         ridge sobre conteos de todos los k-meros con k = 1..5 (1364)
#   kmeros_vent_1a4    ridge sobre conteos de k-meros k = 1..4 en tres ventanas (1020)
#                      (-236 a -61, -60 a -11, -10 a +15 respecto del TSS)
#   cnn                CNN chica sobre one-hot, entrenada desde cero
#                      (multitarea: una salida por replica, perdida enmascarada)
# Embeddings: evaluar_embeddings() recibe {capa: matriz n x d} de un modelo
# y elige la capa DENTRO del entrenamiento de cada fold (error de GCV de
# RidgeCV), nunca mirando el fold de test.
#
# Uso (desde la raiz de TesisDoc):
#   python R/71_exploratorio_embeddings_adn/evaluar.py            # bases sin CNN
#   python R/71_exploratorio_embeddings_adn/evaluar.py --cnn      # bases + CNN
# Salida: figures/exploratorio_embeddings_adn/bases.tsv (+ predicciones_bases.tsv)

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

DATASET = Path("data/processed/embeddings_adn/embeddings_adn_dataset_v1.tsv")
OUT_DIR = Path("figures/exploratorio_embeddings_adn")
REPS = {"Rep 1": "mean_rep1", "Rep 2": "mean_rep2"}
OTRA = {"mean_rep1": "mean_rep2", "mean_rep2": "mean_rep1"}
ADITIVO = ["CGI_frag", "TATA_EPD", "CCAAT_EPD", "YY1_motif", "CRE_motif",
           "other_TSS", "Ebox_CAGGTR", "polyAT_prox"]
ALPHAS = np.logspace(-2, 6, 17)
TSS = 236  # indice 0-based del TSS en la secuencia orientada
VENTANAS = {"distal": (0, TSS - 60), "proximal": (TSS - 60, TSS - 10), "tss": (TSS - 10, 252)}


def cargar(path=DATASET):
    return pd.read_csv(path, sep="\t")


# --- Metricas ------------------------------------------------------------------

def techo(df):
    m = df[["mean_rep1", "mean_rep2"]].dropna()
    return pearsonr(m.mean_rep1, m.mean_rep2)[0]


def metricas(df, pred, col):
    y = df[col].to_numpy()
    ok = ~np.isnan(y) & ~np.isnan(pred)
    r2 = 1 - np.sum((y[ok] - pred[ok]) ** 2) / np.sum((y[ok] - y[ok].mean()) ** 2)
    otra = df[OTRA[col]].to_numpy()
    ok2 = ok & ~np.isnan(otra)
    return {
        "n": int(ok.sum()),
        "pearson": pearsonr(y[ok], pred[ok])[0],
        "spearman": spearmanr(y[ok], pred[ok])[0],
        "r2": r2,
        "pearson_otra_replica": pearsonr(otra[ok2], pred[ok2])[0],
    }


# --- Validacion cruzada por cromosoma --------------------------------------------

def predecir_cv(df, X, col, nuevo_modelo):
    """Predicciones fuera de fold para la replica `col`. X: matriz n x p (o
    dict {capa: matriz} si nuevo_modelo maneja varias capas)."""
    y = df[col].to_numpy()
    pred = np.full(len(df), np.nan)
    info = []
    for f in sorted(df.fold.unique()):
        tr = (df.fold != f).to_numpy() & ~np.isnan(y)
        te = (df.fold == f).to_numpy()
        p, extra = nuevo_modelo(X, tr, te, y)
        pred[te] = p
        info.append({"fold": f, **extra})
    return pred, info


def ajuste_ols(X, tr, te, y):
    m = LinearRegression().fit(X[tr], y[tr])
    return m.predict(X[te]), {}


def ajuste_ridge(X, tr, te, y):
    m = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
    m.fit(X[tr], y[tr])
    return m.predict(X[te]), {"alpha": m[-1].alpha_}


def ajuste_ridge_capas(capas, tr, te, y):
    """Elige la capa con menor error de GCV en el entrenamiento del fold."""
    mejor = None
    for nombre, X in capas.items():
        m = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
        m.fit(X[tr], y[tr])
        score = m[-1].best_score_  # -MSE de leave-one-out (GCV), mayor es mejor
        if mejor is None or score > mejor[0]:
            mejor = (score, nombre, m)
    _, nombre, m = mejor
    return m.predict(capas[nombre][te]), {"capa": nombre, "alpha": m[-1].alpha_}


# --- Caracteristicas de secuencia ----------------------------------------------------

_COD = np.full(256, -1, dtype=np.int64)
for i, b in enumerate("ACGT"):
    _COD[ord(b)] = i


def codificar(seqs):
    return np.stack([_COD[np.frombuffer(s.encode(), dtype=np.uint8)] for s in seqs])


def conteo_kmeros(codigos, k, ini=0, fin=None):
    sub = codigos[:, ini:fin]
    L = sub.shape[1] - k + 1
    idx = np.zeros((sub.shape[0], L), dtype=np.int64)
    for j in range(k):
        idx = idx * 4 + sub[:, j:j + L]
    out = np.zeros((sub.shape[0], 4 ** k), dtype=np.float32)
    filas = np.repeat(np.arange(sub.shape[0]), L)
    np.add.at(out, (filas, idx.ravel()), 1)
    return out


def matriz_kmeros(codigos, ks):
    return np.hstack([conteo_kmeros(codigos, k) for k in ks])


def matriz_kmeros_ventanas(codigos, ks):
    return np.hstack([conteo_kmeros(codigos, k, a, b) for a, b in VENTANAS.values() for k in ks])


# --- CNN desde cero ----------------------------------------------------------------

def cnn_cv(df, epocas_max=40, paciencia=5, semilla=0):
    import torch
    from torch import nn

    torch.manual_seed(semilla)
    np.random.seed(semilla)
    torch.set_num_threads(max(1, torch.get_num_threads()))
    cod = codificar(df.seq)
    onehot = np.eye(4, dtype=np.float32)[cod].transpose(0, 2, 1)  # n x 4 x 252
    Y = df[list(REPS.values())].to_numpy(dtype=np.float32)

    class Red(nn.Module):
        def __init__(self):
            super().__init__()
            def bloque(i, o, k, pool):
                capas = [nn.Conv1d(i, o, k, padding=k // 2), nn.BatchNorm1d(o), nn.ReLU()]
                if pool:
                    capas.append(nn.MaxPool1d(pool))
                return capas
            self.conv = nn.Sequential(*bloque(4, 128, 15, 3), *bloque(128, 128, 7, 3),
                                      *bloque(128, 128, 5, 0), nn.Dropout(0.2))
            self.cabeza = nn.Sequential(nn.Linear(256, 128), nn.ReLU(), nn.Dropout(0.3), nn.Linear(128, 2))

        def forward(self, x):
            h = self.conv(x)
            h = torch.cat([h.mean(-1), h.amax(-1)], dim=1)
            return self.cabeza(h)

    def perdida(pred, y):
        m = ~torch.isnan(y)
        return ((pred[m] - y[m]) ** 2).mean()

    folds = sorted(df.fold.unique())
    pred = np.full_like(Y, np.nan)
    info = []
    for f in folds:
        val_fold = folds[(folds.index(f) + 1) % len(folds)]  # fold interno de validacion
        te = (df.fold == f).to_numpy()
        va = (df.fold == val_fold).to_numpy()
        tr = ~te & ~va
        mu, sd = np.nanmean(Y[tr], 0), np.nanstd(Y[tr], 0)
        Yz = (Y - mu) / sd
        Xt, Yt = torch.tensor(onehot), torch.tensor(Yz)
        red = Red()
        opt = torch.optim.Adam(red.parameters(), lr=1e-3, weight_decay=1e-4)
        mejor, estado, sin_mejora = np.inf, None, 0
        idx_tr = np.where(tr)[0]
        t0 = time.time()
        for ep in range(epocas_max):
            red.train()
            np.random.shuffle(idx_tr)
            for i in range(0, len(idx_tr), 128):
                b = idx_tr[i:i + 128]
                opt.zero_grad()
                l = perdida(red(Xt[b]), Yt[b])
                l.backward()
                opt.step()
            red.eval()
            with torch.no_grad():
                lv = perdida(red(Xt[va]), Yt[va]).item()
            if lv < mejor - 1e-4:
                mejor, sin_mejora = lv, 0
                estado = {k: v.clone() for k, v in red.state_dict().items()}
            else:
                sin_mejora += 1
                if sin_mejora >= paciencia:
                    break
        red.load_state_dict(estado)
        red.eval()
        with torch.no_grad():
            pred[te] = red(Xt[te]).numpy() * sd + mu
        info.append({"fold": f, "epocas": ep + 1, "val_mse": mejor, "seg": round(time.time() - t0)})
        print(f"  CNN fold {f}: {ep + 1} epocas, val MSE {mejor:.3f}, {time.time() - t0:.0f} s", flush=True)
    return {"Rep 1": pred[:, 0], "Rep 2": pred[:, 1]}, info


# --- Orquestacion -------------------------------------------------------------------

def fila(modelo, rep, df, pred, r_techo, extra=None):
    m = metricas(df, pred, REPS[rep])
    return {"modelo": modelo, "replica": rep, **m, "r_entre_replicas": r_techo,
            "frac_techo": m["r2"] / r_techo, **(extra or {})}


def evaluar_bases(df, con_cnn=False):
    r_techo = techo(df)
    cod = codificar(df.seq)
    matrices = {
        "aditivo": (df[ADITIVO].to_numpy(float), ajuste_ols),
        "aditivo_comp": (df[ADITIVO + ["g_c", "cpg_oe"]].to_numpy(float), ajuste_ols),
        "kmeros_1a5": (matriz_kmeros(cod, range(1, 6)), ajuste_ridge),
        "kmeros_vent_1a4": (matriz_kmeros_ventanas(cod, range(1, 5)), ajuste_ridge),
    }
    filas, preds = [], {}
    for nombre, (X, ajuste) in matrices.items():
        t0 = time.time()
        for rep, col in REPS.items():
            pred, info = predecir_cv(df, X, col, ajuste)
            alphas = [i.get("alpha") for i in info if "alpha" in i]
            filas.append(fila(nombre, rep, df, pred, r_techo,
                              {"p": X.shape[1], "alpha_mediana": np.median(alphas) if alphas else np.nan}))
            preds[f"{nombre}|{rep}"] = pred
        print(f"{nombre}: {time.time() - t0:.0f} s", flush=True)
    if con_cnn:
        pcnn, _ = cnn_cv(df)
        for rep in REPS:
            filas.append(fila("cnn", rep, df, pcnn[rep], r_techo))
            preds[f"cnn|{rep}"] = pcnn[rep]
    return pd.DataFrame(filas), pd.DataFrame(preds).assign(seq_id=df.seq_id.values)


def evaluar_embeddings(df, capas, modelo):
    """capas: {nombre_capa: matriz n x d} en el mismo orden de filas que df.
    Devuelve (tabla con la capa elegida dentro de cada fold, tabla por capa)."""
    r_techo = techo(df)
    filas, por_capa = [], []
    for rep, col in REPS.items():
        pred, info = predecir_cv(df, capas, col, ajuste_ridge_capas)
        elegidas = ",".join(str(i["capa"]) for i in info)
        filas.append(fila(modelo, rep, df, pred, r_techo, {"capas_elegidas": elegidas}))
        for nombre, X in capas.items():  # exploratorio: cada capa por separado
            p, _ = predecir_cv(df, X, col, ajuste_ridge)
            por_capa.append(fila(f"{modelo}|{nombre}", rep, df, p, r_techo))
    return pd.DataFrame(filas), pd.DataFrame(por_capa)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cnn", action="store_true", help="incluir la CNN entrenada desde cero")
    ap.add_argument("--dataset", default=str(DATASET))
    ap.add_argument("--out", default=str(OUT_DIR))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df = cargar(a.dataset)
    tabla, preds = evaluar_bases(df, con_cnn=a.cnn)
    tabla.to_csv(out / "bases.tsv", sep="\t", index=False)
    preds.to_csv(out / "predicciones_bases.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(tabla.round(3))


if __name__ == "__main__":
    main()
