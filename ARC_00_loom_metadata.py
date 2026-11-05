# %%
# =============================================================================
# ARC_00_loom_metadata.py   (revision 2)
#
# Read-only. Extracts the cell metadata from the two vendor loom files.
#
# -----------------------------------------------------------------------------
# WHAT REVISION 1 GOT WRONG, AND WHAT ITS OUTPUT ALREADY SETTLED
# -----------------------------------------------------------------------------
#
# THREE BUGS, ALL FIXED HERE
#
#   1. CellID parsed backwards. The real format is BARCODE-1_LIBRARYINDEX, e.g.
#      AAACAGCCACTATGGC-1_1: a 16 bp barcode, the -1 that 10x itself writes, and
#      then _1, the list position Seurat's merge() appends. Revision 1 took the
#      part before the last underscore as the library and the part after as the
#      barcode, which is why barcode lengths came back as {1: 14887}. The
#      LIBRARY INDEX is the trailing integer, and it runs _1.._8 for brain and
#      _1.._12 for colon.
#
#   2. The crash. Embedding, Embeddings_X and Embeddings_Y are h5py COMPOUND
#      dtypes, so dtype.kind is 'V' and revision 1's string test missed them.
#      They entered the frame as numpy void scalars and pandas raised
#      "ufunc 'isnan' not supported" when testing them for NaN on write. They
#      are now expanded into one column per field.
#
#   3. Cluster labels contain LITERAL NEWLINES: "Smooth\nmuscle cells",
#      "Endothelial\nCells-1", "Entero\ncyte-3", "Entero\nchromaffin\ncells",
#      "Transient-\namplifying cells". Unhandled, those corrupt the CSV and
#      break every join. Note the newline sits MID-WORD, so display cleaning
#      gives "Entero cyte-3" while the squash normalizer recovers
#      "enterocyte3" and joins correctly.
#
# WHAT REVISION 1 ALREADY ANSWERED
#
#   THE ARMS. Stimulation against Samples gives
#       SIV/cART/THCCBD   40872, 40874, 40877, 40878   = Arm A
#       SIV/cART/VEH      40873, 40875, 40876, 40879   = Arm B
#   So Arm A is THC/CBD and Arm B is VEH. Figure legends unblocked.
#
#   THE GABAergic CLUSTER. The loom carries 18 brain labels including
#   "GABAergic ADARB2+ Neurons", "GABAergic CCK+ Neurons" and "GABAergic
#   SST+/NPY+ Neurons", and NO glutamatergic cluster. With
#   ATAC_all_clusters.xlsx mapping its CCK+ GABAergic sheet to
#   DA_THCCBDvsVEH_GlutamtergicNeurons.txt, that is two independent sources:
#   the cluster the DA filenames call Glutamtergic is CCK+ GABAergic. Resolved.
#
#   MICROGLIA-2 AND -3. The loom has 18 brain clusters; DA covers 12. The six
#   without DA results are all small, so they were dropped for power.
#
# Jake Lehle, 2026-10
# =============================================================================

# %%
# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

ROOT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome"
LOOMS = {
    "Brain": f"{ROOT}/ATAC_results_Brain/seu_obj_Brain.loom",
    "Colon": f"{ROOT}/ATAC_Results_Colon/seu_obj_Colon.loom",
}

TIDY_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/01_tidy"
OUT_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/02_arc/00_loom_metadata"

MAX_UNIQUE_SHOW = 40
VERBOSE = True


# %%
# -----------------------------------------------------------------------------
# Cell 2: imports
# -----------------------------------------------------------------------------

import os
import re
import sys
import json
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) or ".")
from multiome_common import (
    clean_label, norm_celltype, display_label, parse_cell_id,
    TREATMENT_ARMS, NO_DA_CLUSTERS,
)

try:
    import h5py
except ImportError:
    sys.exit("FATAL: h5py not available. Activate sc_pre_R.")

warnings.filterwarnings("ignore", category=FutureWarning)
os.makedirs(OUT_DIR, exist_ok=True)
NOTES = []

ARM_A = {"40872", "40874", "40877", "40878"}
ARM_B = {"40873", "40875", "40876", "40879"}


def say(m=""):
    if VERBOSE:
        print(m)


def note(m):
    NOTES.append(m)
    say("  [note] " + m)


def decode(arr):
    out = []
    for v in np.asarray(arr).ravel():
        out.append(v.decode("utf-8", errors="replace") if isinstance(v, bytes) else str(v))
    return out


say("=" * 78)
say("ARC_00_loom_metadata rev2")
say("=" * 78)
for t, p in LOOMS.items():
    if not os.path.isfile(p):
        sys.exit(f"FATAL: loom not found for {t}: {p}")
    say(f"  {t:<6} {os.path.basename(p)}  ({os.path.getsize(p)/1e6:.0f} MB)")


# %%
# -----------------------------------------------------------------------------
# Cell 3: read the looms
#
# Compound dtypes are expanded here rather than passed through, which is what
# crashed revision 1. The main UMAP lives in the `Embedding` compound and gets
# named umap_1 / umap_2; Embeddings_X and Embeddings_Y hold alternate
# reductions keyed by name and are kept with a field suffix.
# -----------------------------------------------------------------------------

meta = {}

for tissue, path in LOOMS.items():
    say("")
    say("-" * 78)
    say(tissue)
    say("-" * 78)

    with h5py.File(path, "r") as h5:
        n_cells = h5["matrix"].shape[1]
        say(f"  matrix: {h5['matrix'].shape[0]:,} genes x {n_cells:,} cells")
        if "layers" in h5 and len(h5["layers"]):
            note(f"{tissue}: layers present ({sorted(h5['layers'].keys())}); a "
                 f"peak assay could be hiding there")
        else:
            say("  layers: none, so RNA only and no peak matrix")

        raw = {k: h5["col_attrs"][k][:] for k in sorted(h5["col_attrs"].keys())}
        genes = decode(h5["row_attrs"]["Gene"][:]) if "row_attrs" in h5 and \
            "Gene" in h5["row_attrs"] else []

    df = pd.DataFrame()
    say("")
    say("  column attributes:")
    for k, v in raw.items():
        a = np.asarray(v)
        if a.dtype.names:
            say(f"    {k:<16} COMPOUND with fields {a.dtype.names}")
            for fld in a.dtype.names:
                if k == "Embedding":
                    col = {"_X": "umap_1", "_Y": "umap_2"}.get(fld, f"Embedding_{fld}")
                else:
                    col = f"{k}__{fld}"
                df[col] = np.asarray(a[fld])
            continue
        if a.ndim != 1 or a.shape[0] != n_cells:
            note(f"{tissue}: skipping {k}, shape {a.shape} is not 1-D per cell")
            continue
        say(f"    {k:<16} shape={a.shape} dtype={a.dtype}")
        df[k] = decode(a) if a.dtype.kind in ("S", "O", "U") else a

    df.insert(0, "tissue", tissue)
    df.attrs["genes"] = genes
    meta[tissue] = df
    say("")
    say(f"  tidy frame: {len(df):,} cells x {df.shape[1]} columns")
    say(f"  columns: {list(df.columns)}")
    say("")
    say("  Embedding expanded to umap_1 / umap_2, which is the main embedding")
    say("  and the one to transfer. Embeddings_X and Embeddings_Y are alternate")
    say("  reductions; their '-1' and '1' fields duplicate the main UMAP.")


# %%
# -----------------------------------------------------------------------------
# Cell 4: parse the cell identifiers
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say("Cell identifiers")
say("=" * 78)

for tissue, df in meta.items():
    say("")
    say(f"  {tissue}")
    ids = df["CellID"].astype(str)
    say(f"    examples: {list(ids.head(3))}")

    parsed = ids.map(parse_cell_id)
    df["barcode"] = [p[0] for p in parsed]
    df["tenx_suffix"] = [p[1] for p in parsed]
    df["library_index"] = [p[2] for p in parsed]

    n_bad = df["barcode"].isna().sum()
    if n_bad:
        note(f"{tissue}: {n_bad} identifiers did not parse")

    lens = df["barcode"].dropna().str.len().value_counts()
    say(f"    barcode lengths: {dict(lens)}")
    if set(lens.index) == {16}:
        say("    [ok] all 16 bp, matching the 10x multiome whitelist")
    else:
        note(f"{tissue}: barcodes are not all 16 bp; inspect the parse")

    say(f"    10x suffixes seen: {sorted(df['tenx_suffix'].dropna().unique())}")
    libs = sorted(df["library_index"].dropna().unique())
    say(f"    library indices: {libs}")
    say(f"    {len(libs)} libraries, which should be 8 for brain and 12 for colon")

    dup = df["barcode"].duplicated().sum()
    say(f"    barcodes reused across libraries: {dup:,} of {len(df):,}")
    say("    Reuse is expected: every library draws from the same whitelist, so")
    say("    the library index is what makes an identifier unique.")

    meta[tissue] = df


# %%
# -----------------------------------------------------------------------------
# Cell 5: library index to sample, the aggr ordering
#
# aggr assigns suffixes in CSV order, so this mapping is what our CSV has to
# reproduce.
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say("Library ordering for aggr")
say("=" * 78)

aggr_order = {}

for tissue, df in meta.items():
    say("")
    say(f"  {tissue}")
    sub = df.dropna(subset=["library_index"])
    if not len(sub):
        note(f"{tissue}: no library indices parsed, ordering not recoverable")
        continue

    ct = pd.crosstab(sub["library_index"].astype(int), sub["Samples"])
    say("")
    say("    library index against sample:")
    say(ct.to_string().replace("\n", "\n    "))

    mapping, ambiguous = {}, []
    for idx in ct.index:
        nz = ct.loc[idx][ct.loc[idx] > 0]
        if len(nz) == 1:
            mapping[int(idx)] = nz.index[0]
        else:
            ambiguous.append((int(idx), dict(nz)))

    if ambiguous:
        note(f"{tissue}: {len(ambiguous)} library index/indices span several "
             f"samples, so the ordering is not clean: {ambiguous[:3]}")
    else:
        ordered = [mapping[i] for i in sorted(mapping)]
        aggr_order[tissue] = ordered
        say("")
        say("    RECOVERED. The aggr CSV must list libraries in this order so")
        say("    our barcode suffixes match the vendor's:")
        for i, s in enumerate(ordered, start=1):
            n = int((sub["library_index"].astype(int) == i).sum())
            say(f"      _{i:<3} {s:<10} {n:>7,} cells")


# %%
# -----------------------------------------------------------------------------
# Cell 6: treatment arms
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say("Treatment arms")
say("=" * 78)

arm_resolution = {}

for tissue, df in meta.items():
    say("")
    say(f"  {tissue}")
    say("")
    for v, n in df["Stimulation"].value_counts().items():
        say(f"    {str(v):<24} {n:>8,} cells")

    df["animal"] = df["Samples"].astype(str).str.extract(r"(40\d{3})")[0]
    by_group = df.groupby("Stimulation")["animal"].apply(lambda s: sorted(set(s.dropna())))

    say("")
    say("    animals per group:")
    for g, animals in by_group.items():
        say(f"      {str(g):<24} {animals}")

    say("")
    for g, animals in by_group.items():
        a = set(animals)
        if a == ARM_A:
            arm_resolution[str(g)] = "Arm A"
            say(f"      {g:<24} == Arm A")
        elif a == ARM_B:
            arm_resolution[str(g)] = "Arm B"
            say(f"      {g:<24} == Arm B")
        else:
            say(f"      {g:<24} spans arms: {sorted(a)}")

say("")
say("  RESOLVED:  Arm A is THC/CBD.  Arm B is VEH.")
say("")
say("    SIV/cART/THCCBD and SIV/cART/THC:CBD  (brain and colon spellings)")
say("      40872, 40874, 40877, 40878")
say("    SIV/cART/VEH")
say("      40873, 40875, 40876, 40879")


# %%
# -----------------------------------------------------------------------------
# Cell 7: the pre-infection pairing
#
# The design problem, now checkable against real group labels rather than
# inferred from the sample list.
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say("Pre-infection pairing")
say("=" * 78)

pairing = {}

if "Colon" in meta:
    df = meta["Colon"]
    pre = set(df.loc[df["Stimulation"] == "Pre-infection", "animal"].dropna())
    say("")
    say(f"  animals with a pre-infection colon sample: {sorted(pre)}")
    in_a = sorted(pre & ARM_A)
    in_b = sorted(pre & ARM_B)
    say(f"    in Arm A, THC/CBD:  {in_a}   ({len(in_a)} of {len(pre)})")
    say(f"    in Arm B, VEH:      {in_b}   ({len(in_b)} of {len(pre)})")
    pairing = {"pre_animals": sorted(pre), "pre_in_THCCBD": in_a, "pre_in_VEH": in_b}

    say("")
    say("  Both versus-pre contrasts draw on the SAME four-animal baseline, but")
    say(f"  {len(in_a)} of those four are THC/CBD and {len(in_b)} is VEH. So")
    say("  THCCBD_vs_Pre is largely a within-animal paired comparison while")
    say("  VEH_vs_Pre is almost entirely unpaired. The two are not comparable")
    say("  to each other, and differences between them partly reflect pairing")
    say("  structure rather than treatment.")
    say("")
    say("  VEH_vs_THCCBD remains the only colon contrast with a balanced")
    say("  4-versus-4 design and should lead.")

    say("")
    say("  cells per group:")
    for g, n in df["Stimulation"].value_counts().items():
        say(f"    {str(g):<24} {n:>8,}")
    say("")
    say("  Pre-infection carries roughly twice the cells of either treatment")
    say("  arm. That inflates a cell-level test asymmetrically and gives uneven")
    say("  library sizes for the pseudobulk recompute, so it belongs in the")
    say("  methods either way.")


# %%
# -----------------------------------------------------------------------------
# Cell 8: cluster labels, cleaned and normalized
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say("Cluster labels")
say("=" * 78)

for tissue, df in meta.items():
    df["Clusters_raw"] = df["Clusters"]
    df["Clusters_display"] = df["Clusters"].map(clean_label)
    df["celltype"] = df["Clusters"].map(norm_celltype)
    df["celltype_display"] = df["celltype"].map(display_label)
    meta[tissue] = df

    n_nl = df["Clusters_raw"].astype(str).str.contains("\n").sum()
    say("")
    say(f"  {tissue}: {df['Clusters_raw'].nunique()} labels, "
        f"{n_nl:,} cells carry a label containing a literal newline")

    tab = (df.groupby(["Clusters_display", "celltype", "celltype_display"])
             .size().reset_index(name="n").sort_values("n", ascending=False))
    say("")
    say(f"    {'loom label':<30} {'token':<26} {'n':>7}")
    say("    " + "-" * 66)
    for r in tab.itertuples():
        flag = "  <-- UNRECOGNIZED" if str(r.celltype).startswith("UNKNOWN:") else ""
        say(f"    {r.Clusters_display:<30} {r.celltype:<26} {r.n:>7,}{flag}")

    unk = tab[tab["celltype"].astype(str).str.startswith("UNKNOWN:")]
    for r in unk.itertuples():
        note(f"{tissue}: cluster '{r.Clusters_display}' does not normalize; add "
             f"it to CELLTYPE_ALIASES in multiome_common.py")


# %%
# -----------------------------------------------------------------------------
# Cell 9: loom clusters against the DA results
#
# Which annotated cell types have no differential accessibility, and therefore
# cannot appear in any DA figure.
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say("Loom clusters against the DA results")
say("=" * 78)

da_tokens = {}
peaks_path = os.path.join(TIDY_DIR, "peaks_da.csv")
if os.path.isfile(peaks_path):
    try:
        pk = pd.read_csv(peaks_path, usecols=["tissue", "celltype"])
        for t, g in pk.groupby("tissue"):
            da_tokens[t] = set(g["celltype"].unique())
        say(f"  read DA cell types from {os.path.basename(peaks_path)}")
    except Exception as exc:
        note(f"could not read {peaks_path} ({exc})")
else:
    note(f"{peaks_path} not found; run MULTIOME_01_build_tidy.py to enable "
         f"this comparison")

no_da = {}
for tissue, df in meta.items():
    if tissue not in da_tokens:
        continue
    loom_set = set(df["celltype"].unique())
    da_set = da_tokens[tissue]
    missing = sorted(loom_set - da_set)
    extra = sorted(da_set - loom_set)

    say("")
    say(f"  {tissue}: {len(loom_set)} loom clusters, {len(da_set)} with DA results")

    if missing:
        say("")
        say("    annotated but NO DA results, so they cannot appear in a DA figure:")
        for t in missing:
            n = int((df["celltype"] == t).sum())
            say(f"      {display_label(t):<30} {n:>7,} cells")
        no_da[tissue] = {t: int((df["celltype"] == t).sum()) for t in missing}

        big = {t: n for t, n in no_da[tissue].items() if n >= 300}
        if big:
            note(f"{tissue}: these have no DA results despite a substantial cell "
                 f"count, which is worth asking about rather than assuming was "
                 f"a power decision: "
                 + ", ".join(f"{display_label(t)} ({n:,})" for t, n in big.items()))

    if extra:
        note(f"{tissue}: DA results exist for tokens absent from the loom: {extra}")

say("")
say("  Recorded expectation in multiome_common.NO_DA_CLUSTERS:")
for t, d in NO_DA_CLUSTERS.items():
    say(f"    {t}: {d}")


# %%
# -----------------------------------------------------------------------------
# Cell 10: write
# -----------------------------------------------------------------------------

say("")
say("=" * 78)
say(f"Writing to {OUT_DIR}")
say("=" * 78)

KEEP = ["tissue", "CellID", "barcode", "tenx_suffix", "library_index",
        "Samples", "animal", "Stimulation",
        "Clusters_raw", "Clusters_display", "celltype", "celltype_display",
        "umap_1", "umap_2", "nGene", "nUMI"]

for tissue, df in meta.items():
    cols = [c for c in KEEP if c in df.columns]
    out = df[cols].copy()
    # Clusters_raw holds embedded newlines; quoting keeps the CSV valid but the
    # display and token columns are what downstream code should use.
    p = os.path.join(OUT_DIR, f"cell_metadata_{tissue}.csv")
    out.to_csv(p, index=False, lineterminator="\n")
    say(f"  wrote cell_metadata_{tissue}.csv  ({len(out):,} cells, {len(cols)} cols)")

combined = pd.concat(
    [df[[c for c in KEEP if c in df.columns]] for df in meta.values()],
    ignore_index=True)
combined.to_csv(os.path.join(OUT_DIR, "cell_metadata_all.csv"),
                index=False, lineterminator="\n")
say(f"  wrote cell_metadata_all.csv  ({len(combined):,} cells)")

# The barcode-to-annotation map the rebuild consumes, keyed so it can be joined
# after aggr.
for tissue, df in meta.items():
    cols = ["library_index", "barcode", "Samples", "Stimulation",
            "celltype", "celltype_display", "umap_1", "umap_2"]
    cols = [c for c in cols if c in df.columns]
    p = os.path.join(OUT_DIR, f"annotation_transfer_{tissue}.csv")
    df[cols].to_csv(p, index=False, lineterminator="\n")
    say(f"  wrote annotation_transfer_{tissue}.csv")

payload = {
    "generated": pd.Timestamp.now().isoformat(),
    "looms": LOOMS,
    "n_cells": {t: int(len(d)) for t, d in meta.items()},
    "cell_id_format": "BARCODE-<tenx>_<library_index>, Seurat merge() convention",
    "aggr_library_order": aggr_order,
    "arm_resolution": arm_resolution,
    "arms": {"Arm A": "THC/CBD", "Arm B": "VEH"},
    "arm_membership": {"Arm A": sorted(ARM_A), "Arm B": sorted(ARM_B)},
    "pre_infection_pairing": pairing,
    "clusters": {
        t: sorted(set(d["celltype_display"])) for t, d in meta.items()
    },
    "clusters_without_da": no_da,
    "gabaergic_resolution": (
        "The cluster every DA filename calls GlutamtergicNeurons is CCK+ "
        "GABAergic. Two independent sources agree: the loom carries GABAergic "
        "ADARB2+, CCK+ and SST+/NPY+ Neurons with no glutamatergic cluster at "
        "all, and ATAC_all_clusters.xlsx maps its CCK+ GABAergic sheet to that "
        "DA file. The join token stays Glutamatergic_FILENAME; figures use "
        "DISPLAY_LABELS."
    ),
    "notes": NOTES,
}
with open(os.path.join(OUT_DIR, "loom_metadata_summary.json"), "w") as fh:
    json.dump(payload, fh, indent=2, default=str)
say("  wrote loom_metadata_summary.json")

say("")
say(f"  {len(NOTES)} note(s)")
for n in NOTES:
    say(f"    - {n}")

say("")
say("Resolved by this run:")
say("  Arm A is THC/CBD, Arm B is VEH")
say("  The DA token GlutamtergicNeurons is CCK+ GABAergic")
say("  aggr library ordering recovered from the Seurat merge() indices")
say("")
say("Next:")
say("  python ARC_01_make_libraries.py")
say("")
say("Done. Nothing was modified.")
