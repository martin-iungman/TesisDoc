#!/usr/bin/env python3
"""Filtra hg38.fa.align.gz (RepeatMasker, UCSC) a los TE que se superponen
con la library de promotores y calcula la divergencia de Kimura (K2P) de
cada copia contra su consenso, con y sin la correccion CpG de
calcDivergenceFromAlign.pl:
  - sitio CpG del consenso (CG) con UNA transicion en la copia (TG o CA):
    cuenta 1/10 de una transicion;
  - sitio CpG con las dos transiciones (TA): cuenta como UNA transicion;
  - transversiones: cuentan normal.
K2P = -0.5*ln(1-2p-q) - 0.25*ln(1-2q), p = transiciones/L, q =
transversiones/L, L = columnas sin gaps (solo A/C/G/T). Se devuelve en
% (comparable con perc_div del .out).

El alineamiento se lee en streaming (el .gz de 2.5GB no se descomprime a
disco). Salida: una fila por (copia de TE, promotor) con superposicion
> 10 pb, excluyendo Simple_repeat/Low_complexity (misma regla que
build_prom_features.R).

Uso (desde la raiz de TesisDoc):
  python R/00_prom_features/filter_rmsk_align.py <align.gz> data/raw/library.bed data/processed/te_hits_kimura.tsv
"""
import bisect
import gzip
import math
import sys
from collections import defaultdict

align_path, bed_path, out_path = sys.argv[1:4]
TRANS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
LCR = ("Simple_repeat", "Low_complexity")

lib = defaultdict(list)
with open(bed_path) as f:
    for line in f:
        t = line.split("\t")
        if len(t) >= 4:
            lib[t[0]].append((int(t[1]) + 1, int(t[2]), t[3].strip()))  # 1-based inclusive
max_w = 0
for chrom in lib:
    lib[chrom].sort()
    max_w = max(max_w, max(e - s + 1 for s, e, _ in lib[chrom]))
starts = {c: [s for s, _, _ in v] for c, v in lib.items()}


def overlaps(chrom, hs, he):
    if chrom not in lib:
        return []
    i = bisect.bisect_left(starts[chrom], hs - max_w)
    out = []
    for s, e, sid in lib[chrom][i:]:
        if s > he:
            break
        ov = min(e, he) - max(s, hs) + 1
        if ov > 10:
            out.append((sid, ov))
    return out


def kimura(g, c):
    """Devuelve (K_cpg, K_nocpg, L, n_cpg_sites) en %."""
    L = 0
    tv = 0
    ti_nocpg = 0
    ti_cpg = 0.0
    n_cpg = 0
    n = len(g)
    k = 0
    while k < n:
        a, b = g[k], c[k]
        if a == "-" or b == "-" or a not in "ACGT" or b not in "ACGT":
            k += 1
            continue
        if b == "C" and k + 1 < n and c[k + 1] == "G" and g[k + 1] in "ACGT":
            a2, b2 = g[k + 1], "G"
            L += 2
            n_cpg += 1
            t1 = a != b and (b, a) in TRANS
            t2 = a2 != b2 and (b2, a2) in TRANS
            v1 = a != b and not t1
            v2 = a2 != b2 and not t2
            tv += v1 + v2
            ti_nocpg += t1 + t2
            if t1 and t2:
                ti_cpg += 1.0
            elif t1 or t2:
                ti_cpg += 0.1
            k += 2
            continue
        L += 1
        if a != b:
            if (b, a) in TRANS:
                ti_nocpg += 1
                ti_cpg += 1.0
            else:
                tv += 1
        k += 1
    if L == 0:
        return None, None, 0, n_cpg

    def k2p(ti):
        p, q = ti / L, tv / L
        x, y = 1 - 2 * p - q, 1 - 2 * q
        if x <= 0 or y <= 0:
            return None
        return 100 * (-0.5 * math.log(x) - 0.25 * math.log(y))

    return k2p(ti_cpg), k2p(ti_nocpg), L, n_cpg


def fmt(x):
    return "NA" if x is None else ("%.3f" % x)


out = open(out_path, "w")
out.write("seq_id\tchr\thit_start\thit_end\toverlap\tname\tclass\tperc_div_hdr\tkimura_cpg\tkimura_nocpg\taln_len\tn_cpg\n")

keep = False
cur = None
G, C = [], []
n_hits = n_kept = 0
last_chr = None


def flush():
    global n_kept
    if cur is None:
        return
    g, c = "".join(G).upper(), "".join(C).upper()
    if len(g) != len(c):  # bloque truncado/mal formado: se descarta
        print("bloque con largos distintos, descartado: %s %s %s" % (cur["chr"], cur["s"], cur["name"]), file=sys.stderr)
        return
    kc, kn, L, ncpg = kimura(g, c)
    for sid, ov in cur["ovs"]:
        out.write("\t".join(map(str, [sid, cur["chr"], cur["s"], cur["e"], ov, cur["name"], cur["cls"], cur["div"], fmt(kc), fmt(kn), L, ncpg])) + "\n")
    n_kept += 1


with gzip.open(align_path, "rb") as f:
    for raw in f:
        if raw[:1].isdigit():
            if keep:
                flush()
            keep = False
            cur = None
            G, C = [], []
            n_hits += 1
            t = raw.decode("latin-1").split()
            chrom, hs, he = t[4], int(t[5]), int(t[6])
            if chrom != last_chr:
                print("%s (hits leidos: %d, TE guardados: %d)" % (chrom, n_hits, n_kept), file=sys.stderr, flush=True)
                last_chr = chrom
            tok = t[9] if t[8] == "C" else t[8]
            name, _, cls = tok.partition("#")
            if cls.startswith(LCR) or cls in LCR:
                continue
            ovs = overlaps(chrom, hs, he)
            if ovs:
                keep = True
                cur = {"chr": chrom, "s": hs, "e": he, "name": name, "cls": cls, "div": t[1], "ovs": ovs}
            continue
        if not keep:
            continue
        line = raw.decode("latin-1")
        if not line.strip() or line.startswith(("Matrix", "Transitions", "Gap_init", "Kimura")):
            continue
        t = line.split()
        # genoma: "  chrN  start SEQ end"; consenso: "[C] nombre  start SEQ end"
        # (en hebra + la linea del consenso TAMBIEN empieza con espacios, asi
        # que se distingue por estructura y no por indentacion)
        if len(t) == 4 and t[0] == cur["chr"]:
            G.append(t[2])
        elif len(t) >= 4 and t[-1].isdigit() and t[-3].isdigit():
            C.append(t[-2])
if keep:
    flush()
out.close()
print("listo: %d hits leidos, %d TE guardados" % (n_hits, n_kept), file=sys.stderr)
