# EXPLORATORIO - embeddings de ADN: mejoras del modelo (corre local, sin GPU)
#
# Parte de los embeddings ya extraidos en Colab (no los recalcula) y prueba,
# todo en escala de PERCENTIL dentro de cada replica (techo = Pearson entre
# los percentiles de las dos replicas, ~0.69; ver handoff.md, "Escala"):
#
#   1. Bases en percentiles, incluida una base mas exigente (k-meros por
#      ventana con k <= 5).
#   2. Embeddings con la capa que gano en la escala cruda (elegida dentro del
#      entrenamiento de cada fold en todos los folds, ver embeddings.tsv):
#      HyenaDNA capa 2, NT v2 capa 11, DNABERT-2 capas 3 (ventanas) y 6 (media).
#   3. Combinaciones: dos capas de DNABERT-2; DNABERT-2 + k-meros + aditivo;
#      los tres modelos de lenguaje juntos.
#   4. Cabezal no lineal: gradient boosting sobre DNABERT-2 reducido por PCA
#      (ajustado dentro de cada fold) + caracteristicas del aditivo.
#   5. IC 95% por bootstrap de la diferencia de R2 de cada modelo contra los
#      k-meros por ventana y contra el aditivo, remuestreando promotores sobre
#      las predicciones fuera de fold (mismos promotores para los dos modelos).
#
# Mismo split por cromosoma y mismo ridge que evaluar.py. Cada modelo guarda
# sus predicciones fuera de fold en mejoras/predicciones/<modelo>.tsv; si ese
# archivo existe, el modelo no se vuelve a correr (usar --rehacer para forzar).
#
# Requiere: data/processed/embeddings_adn/embeddings_adn_dataset_v1.tsv (preparar_dataset.py)
# y data/processed/embeddings_adn/emb_{hyenadna,nt_v2_100m,dnabert2}.npz
# (bajados de la carpeta de Drive donde los dejo el notebook de Colab).
# Paquetes: numpy, pandas, scipy, scikit-learn.
#
# Uso, desde la raiz de TesisDoc:
#   python R/71_exploratorio_embeddings_adn/mejoras.py
#   python R/71_exploratorio_embeddings_adn/mejoras.py --solo dnabert2_vent_c3,kmeros_vent_1a5
#   python R/71_exploratorio_embeddings_adn/mejoras.py --listar
#   python R/71_exploratorio_embeddings_adn/mejoras.py --objetivo promedio
# Salida: figures/exploratorio_embeddings_adn/mejoras/{resultados,bootstrap}.tsv
#
# --objetivo promedio: solo los promotores medidos en las dos replicas (8889).
# Percentil dentro de cada replica (calculado en ese subconjunto) y se entrena
# y evalua contra el PROMEDIO de los dos percentiles. Techo = confiabilidad del
# promedio por Spearman-Brown, 2r / (1 + r), con r = Pearson entre los
# percentiles de las dos replicas. Salida en mejoras_promedio/ (mismo formato).

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluar as ev  # noqa: E402

EMB_DIR = Path("data/processed/embeddings_adn")
OUT = ev.OUT_DIR / "mejoras"
OUT_PROM = ev.OUT_DIR / "mejoras_promedio"
N_BOOT = 1000
SEMILLA = 0


# --- Datos ---------------------------------------------------------------------

def dataset_percentil():
    df = ev.cargar()
    for col in ev.REPS.values():
        df[col] = 100 * df[col].rank(pct=True)  # NaN se mantiene
    return df


def dataset_promedio():
    """Promotores medidos en las dos replicas; columna pct_prom = promedio de
    los percentiles de cada replica (calculados en este subconjunto). Devuelve
    tambien la mascara sobre el dataset completo (para recortar embeddings)."""
    df = ev.cargar()
    filas = (df.mean_rep1.notna() & df.mean_rep2.notna()).to_numpy()
    sub = df[filas].reset_index(drop=True)
    for col in ev.REPS.values():
        sub[f"pct_{col}"] = 100 * sub[col].rank(pct=True)
    sub["pct_prom"] = (sub.pct_mean_rep1 + sub.pct_mean_rep2) / 2
    return sub, filas


_cache_npz = {}


def emb(modelo, clave):
    """Matriz float32 de un embedding guardado (clave: 'capa3_ventanas', etc.)."""
    path = EMB_DIR / f"emb_{modelo}.npz"
    if not path.exists():
        sys.exit(f"Falta {path}: bajarlo de la carpeta de Drive (ver handoff.md)")
    if path not in _cache_npz:
        _cache_npz[path] = np.load(path)
    return _cache_npz[path][clave].astype(np.float32)


# --- Ajustes -------------------------------------------------------------------

def ajuste_hgb_pca(X, tr, te, y, n_pca=128, n_pca_cols=None):
    """Gradient boosting. Si n_pca_cols esta dado, las primeras n_pca_cols
    columnas (el embedding) se reducen por PCA dentro del fold y el resto
    (caracteristicas del aditivo) entra tal cual."""
    Xe, Xr = X[:, :n_pca_cols], X[:, n_pca_cols:]
    red = make_pipeline(StandardScaler(), PCA(n_components=min(n_pca, Xe.shape[1]), random_state=SEMILLA))
    Ztr = np.hstack([red.fit_transform(Xe[tr]), Xr[tr]])
    Zte = np.hstack([red.transform(Xe[te]), Xr[te]])
    m = HistGradientBoostingRegressor(learning_rate=0.05, max_iter=1000, max_leaf_nodes=31,
                                      l2_regularization=1.0, early_stopping=True,
                                      validation_fraction=0.15, n_iter_no_change=30,
                                      random_state=SEMILLA)
    m.fit(Ztr, y[tr])
    return m.predict(Zte), {"iteraciones": m.n_iter_}


# --- Catalogo de modelos ---------------------------------------------------------
# Cada entrada: (descripcion, funcion que arma X a partir de df, funcion de ajuste)

def _aditivo(df):
    return df[ev.ADITIVO].to_numpy(float)


def _aditivo_comp(df):
    return df[ev.ADITIVO + ["g_c", "cpg_oe"]].to_numpy(float)


def _kmeros_vent(ks):
    return lambda df: ev.matriz_kmeros_ventanas(ev.codificar(df.seq), ks)


def _concat(*fs):
    return lambda df: np.hstack([f(df) for f in fs])


def _emb(modelo, clave):
    return lambda df: emb(modelo, clave)


def _hgb_dnabert(df):
    return np.hstack([emb("dnabert2", "capa3_ventanas"), _aditivo_comp(df)])


MODELOS = {
    # 1. bases
    "aditivo": ("8 elementos del aditivo (OLS)", _aditivo, ev.ajuste_ols),
    "aditivo_comp": ("aditivo + G+C y CpG o/e (OLS)", _aditivo_comp, ev.ajuste_ols),
    "kmeros_vent_1a4": ("k-meros k=1..4 en tres ventanas", _kmeros_vent(range(1, 5)), ev.ajuste_ridge),
    "kmeros_vent_1a5": ("k-meros k=1..5 en tres ventanas", _kmeros_vent(range(1, 6)), ev.ajuste_ridge),
    # 2. embeddings, capa fija
    "hyenadna_vent_c2": ("HyenaDNA, ventanas, capa 2", _emb("hyenadna", "capa2_ventanas"), ev.ajuste_ridge),
    "nt_vent_c11": ("NT v2 100M, ventanas, capa 11", _emb("nt_v2_100m", "capa11_ventanas"), ev.ajuste_ridge),
    "dnabert2_media_c6": ("DNABERT-2, media, capa 6", _emb("dnabert2", "capa6_media"), ev.ajuste_ridge),
    "dnabert2_vent_c3": ("DNABERT-2, ventanas, capa 3", _emb("dnabert2", "capa3_ventanas"), ev.ajuste_ridge),
    # 3. combinaciones
    "dnabert2_vent_c3c6": ("DNABERT-2, ventanas, capas 3 + 6",
                           _concat(_emb("dnabert2", "capa3_ventanas"), _emb("dnabert2", "capa6_ventanas")),
                           ev.ajuste_ridge),
    "dnabert2_kmeros_aditivo": ("DNABERT-2 ventanas c3 + k-meros 1..4 por ventana + aditivo_comp",
                                _concat(_emb("dnabert2", "capa3_ventanas"), _kmeros_vent(range(1, 5)), _aditivo_comp),
                                ev.ajuste_ridge),
    "tres_modelos_vent": ("HyenaDNA c2 + NT c11 + DNABERT-2 c3, ventanas",
                          _concat(_emb("hyenadna", "capa2_ventanas"), _emb("nt_v2_100m", "capa11_ventanas"),
                                  _emb("dnabert2", "capa3_ventanas")),
                          ev.ajuste_ridge),
    # 4. no lineal
    "hgb_dnabert2_aditivo": ("gradient boosting: PCA(128) de DNABERT-2 ventanas c3 + aditivo_comp",
                             _hgb_dnabert,
                             lambda X, tr, te, y: ajuste_hgb_pca(
                                 X, tr, te, y, n_pca_cols=emb("dnabert2", "capa3_ventanas").shape[1])),
}

# Modelos entrenados fuera de este script (finetuning.py en Colab): solo se leen
# sus predicciones si estan en predicciones/ y entran al resumen y al bootstrap.
EXTERNOS = {
    "cnn_ajustada": "CNN con hiperparametros elegidos en validacion interna, promedio de 3 semillas (finetuning.py)",
    "cnn_ajustada_s0": "CNN ajustada, una semilla (finetuning.py)",
    "dnabert2_ft": "DNABERT-2 con fine-tuning, pooling por ventanas, promedio de 3 semillas (finetuning.py)",
    "dnabert2_ft_s0": "DNABERT-2 con fine-tuning, una semilla (finetuning.py)",
}

# cnn_ajustada y dnabert2_vent_c3c6: para comparar el fine-tuning contra la base
# supervisada y contra el mejor DNABERT-2 congelado (se usan solo si existen).
REFERENCIAS = ["kmeros_vent_1a4", "aditivo", "cnn_ajustada", "dnabert2_vent_c3c6"]


# --- Correr -----------------------------------------------------------------------

def correr_modelo(nombre, df, objetivos, pred_dir, df_completo=None, filas=None, rehacer=False):
    """objetivos: {etiqueta: columna de df}. Si filas esta dado, X se arma sobre
    df_completo (mismo orden que los embeddings) y se recorta a esas filas."""
    path = pred_dir / f"{nombre}.tsv"
    if path.exists() and not rehacer:
        print(f"{nombre}: ya estaba, se usa {path}")
        return pd.read_csv(path, sep="\t")
    desc, armar_X, ajuste = MODELOS[nombre]
    t0 = time.time()
    X = armar_X(df if filas is None else df_completo)
    if filas is not None:
        X = X[filas]
    out = pd.DataFrame({"seq_id": df.seq_id})
    for rep, col in objetivos.items():
        pred, _ = ev.predecir_cv(df, X, col, ajuste)
        out[rep] = pred
    out.to_csv(path, sep="\t", index=False)
    print(f"{nombre}: {X.shape[1]} variables, {time.time() - t0:.0f} s", flush=True)
    return out


def r2(y, p):
    return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)


def fila_promedio(modelo, df, pred, r_techo, extra=None):
    y = df.pct_prom.to_numpy()
    return {"modelo": modelo, "replica": "Promedio", "n": len(y),
            "pearson": pearsonr(y, pred)[0], "spearman": spearmanr(y, pred)[0], "r2": r2(y, pred),
            "pearson_pct_rep1": pearsonr(df.pct_mean_rep1, pred)[0],
            "pearson_pct_rep2": pearsonr(df.pct_mean_rep2, pred)[0],
            "techo_spearman_brown": r_techo, "frac_techo": r2(y, pred) / r_techo, **(extra or {})}


def bootstrap(df, preds, r_techo, objetivos):
    rng = np.random.default_rng(SEMILLA)
    filas = []
    for rep, col in objetivos.items():
        ok = df[col].notna().to_numpy()
        y = df[col].to_numpy()[ok]
        idx = rng.integers(0, len(y), size=(N_BOOT, len(y)))
        for ref in REFERENCIAS:
            if ref not in preds:
                continue
            pr = preds[ref][rep].to_numpy()[ok]
            for nombre, p in preds.items():
                if nombre == ref:
                    continue
                pm = p[rep].to_numpy()[ok]
                dif = np.array([r2(y[i], pm[i]) - r2(y[i], pr[i]) for i in idx])
                obs = r2(y, pm) - r2(y, pr)
                filas.append({"modelo": nombre, "referencia": ref, "replica": rep,
                              "dif_r2": obs, "ic_inf": np.quantile(dif, 0.025),
                              "ic_sup": np.quantile(dif, 0.975),
                              "dif_frac_techo": obs / r_techo})
    return pd.DataFrame(filas)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo", help="modelos separados por coma (default: todos)")
    ap.add_argument("--rehacer", action="store_true", help="recalcular aunque existan predicciones")
    ap.add_argument("--listar", action="store_true", help="mostrar el catalogo y salir")
    ap.add_argument("--objetivo", choices=["replica", "promedio"], default="replica",
                    help="replica: cada replica por separado (default); "
                         "promedio: promedio de percentiles, promotores medidos en ambas")
    a = ap.parse_args()
    if a.listar:
        for k, (d, _, _) in MODELOS.items():
            print(f"{k:28s} {d}")
        return
    nombres = a.solo.split(",") if a.solo else list(MODELOS)
    if a.objetivo == "replica":
        out_dir, df_completo, filas_X = OUT, None, None
        df = dataset_percentil()
        objetivos = ev.REPS
        r_techo = ev.techo(df)
        print(f"Techo en percentiles (Pearson entre replicas): {r_techo:.3f}")
    else:
        out_dir, df_completo = OUT_PROM, ev.cargar()
        df, filas_X = dataset_promedio()
        objetivos = {"Promedio": "pct_prom"}
        r = pearsonr(df.pct_mean_rep1, df.pct_mean_rep2)[0]
        r_techo = 2 * r / (1 + r)
        print(f"n = {len(df)}; r entre percentiles = {r:.3f}; techo Spearman-Brown = {r_techo:.3f}")
    pred_dir = out_dir / "predicciones"
    pred_dir.mkdir(parents=True, exist_ok=True)

    for n in nombres:
        if n in EXTERNOS:
            continue  # se entrena en finetuning.py; aca solo se leen sus predicciones
        correr_modelo(n, df, objetivos, pred_dir, df_completo, filas_X, a.rehacer)

    # Resumen con TODAS las predicciones disponibles (tambien las de corridas previas)
    desc = {**{k: v[0] for k, v in MODELOS.items()}, **EXTERNOS}
    preds = {p.stem: pd.read_csv(p, sep="\t") for p in sorted(pred_dir.glob("*.tsv")) if p.stem in desc}
    filas = []
    for n, p in preds.items():
        assert (p.seq_id.values == df.seq_id.values).all(), f"{n}: seq_id no coincide con el dataset"
        extra = {"descripcion": desc[n]}
        if a.objetivo == "replica":
            for rep in ev.REPS:
                filas.append(ev.fila(n, rep, df, p[rep].to_numpy(), r_techo, extra))
        else:
            filas.append(fila_promedio(n, df, p["Promedio"].to_numpy(), r_techo, extra))
    res = pd.DataFrame(filas)
    res.to_csv(out_dir / "resultados.tsv", sep="\t", index=False)
    boot = bootstrap(df, preds, r_techo, objetivos)
    boot.to_csv(out_dir / "bootstrap.tsv", sep="\t", index=False)

    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(res.pivot_table(index="modelo", columns="replica", values=["r2", "frac_techo"])
              .round(3).sort_values(("r2", next(iter(objetivos)))))
        print(boot[boot.referencia == "kmeros_vent_1a4"].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
