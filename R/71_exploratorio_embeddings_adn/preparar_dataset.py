# EXPLORATORIO - embeddings de modelos de lenguaje de ADN (paso 0: dataset)
#
# Arma una tabla unica, una fila por promotor (type == "promoter"), con:
#   - la secuencia ensayada de 252 pb, orientada en el sentido de la
#     transcripcion (TSS = posicion 237, 1-based), tomada de
#     transcriptional_library/Library_data/data/promoters_wo_dupl.fa;
#   - la actividad media de cada replica (activity_stats_highconf.tsv,
#     columna mean; NA si el promotor no paso el filtro en esa replica);
#   - las caracteristicas del modelo aditivo de R/38_aditividad (cinco
#     positivas y tres negativas) mas G+C y CpG o/e continuos;
#   - un fold de validacion cruzada (0-4) asignado por CROMOSOMA.
#
# Por que por cromosoma: promotores alternativos de un mismo gen, familias
# de transposones y duplicaciones comparten secuencia. Si un pariente queda
# en entrenamiento y el otro en test, el modelo acierta por memoria. Todos
# los promotores de un gen estan en el mismo cromosoma, asi que el split
# por cromosoma tambien agrupa por gen. Los cromosomas se reparten en 5
# folds de tamano parecido (greedy, del mas grande al mas chico); chrY (1
# promotor) va con chrX. La asignacion queda fija en el dataset, para que
# todos los modelos (bases y embeddings) usen exactamente el mismo split.
#
# Salida: data/processed/embeddings_adn/embeddings_adn_dataset_v1.tsv
# Correr desde la raiz de TesisDoc:  python R/71_exploratorio_embeddings_adn/preparar_dataset.py

import sys
from pathlib import Path

import pandas as pd

FASTA = Path("../transcriptional_library/Library_data/data/promoters_wo_dupl.fa")
PROM_DF = Path("data/processed/prom_df.tsv")
ACT = Path("data/processed/activity_stats_highconf.tsv")
OUT = Path("data/processed/embeddings_adn/embeddings_adn_dataset_v1.tsv")
N_FOLDS = 5

POSITIVOS = ["CGI_frag", "TATA_EPD", "CCAAT_EPD", "YY1_motif", "CRE_motif"]
NEGATIVOS = ["other_TSS", "Ebox_CAGGTR", "polyAT_prox"]
CONTINUAS = ["g_c", "cpg_oe"]


def leer_fasta(path):
    seqs, k = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                k = line[1:]
            elif k is not None:
                seqs[k] = seqs.get(k, "") + line.upper()
    return seqs


def asignar_folds(n_por_crom, n_folds):
    carga = [0] * n_folds
    fold = {}
    for crom, n in sorted(n_por_crom.items(), key=lambda x: -x[1]):
        if crom == "chrY":
            continue
        i = carga.index(min(carga))
        fold[crom] = i
        carga[i] += n
    if "chrY" in n_por_crom:
        fold["chrY"] = fold["chrX"]
    return fold


def main():
    seqs = leer_fasta(FASTA)
    prom = pd.read_csv(PROM_DF, sep="\t", low_memory=False)
    prom = prom.loc[prom["type"] == "promoter",
                    ["seq_id", "gene_sym", "chromosome"] + POSITIVOS + NEGATIVOS + CONTINUAS].copy()
    faltan = set(prom.seq_id) - set(seqs)
    if faltan:
        sys.exit(f"{len(faltan)} promotores sin secuencia en el FASTA")
    prom["seq"] = prom.seq_id.map(seqs)
    assert prom.seq.str.len().eq(252).all() and prom.seq.str.fullmatch("[ACGT]+").all()

    act = pd.read_csv(ACT, sep="\t")
    act = act.pivot_table(index="seq_id", columns="rep", values="mean")
    act.columns = ["mean_rep1", "mean_rep2"]
    prom = prom.merge(act, left_on="seq_id", right_index=True, how="left")

    folds = asignar_folds(prom.chromosome.value_counts().to_dict(), N_FOLDS)
    prom["fold"] = prom.chromosome.map(folds)
    for c in POSITIVOS + NEGATIVOS:
        prom[c] = prom[c].astype(int)

    # Solo promotores medidos en al menos una replica
    prom = prom.dropna(subset=["mean_rep1", "mean_rep2"], how="all").reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    prom.to_csv(OUT, sep="\t", index=False)

    print(f"{len(prom)} promotores -> {OUT}")
    print("Medidos: Rep 1", prom.mean_rep1.notna().sum(), "| Rep 2", prom.mean_rep2.notna().sum(),
          "| ambas", prom[["mean_rep1", "mean_rep2"]].notna().all(axis=1).sum())
    print(prom.groupby("fold").agg(n=("seq_id", "size"),
                                   cromosomas=("chromosome", lambda x: ",".join(sorted(set(x))))))


if __name__ == "__main__":
    main()
