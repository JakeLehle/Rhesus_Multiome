# %%
# =============================================================================
# MULTIOME_00_inventory.py
#
# Read-only inventory and schema diagnostic for the THC/CBD rhesus multiome
# collaboration (Genohub-processed brain and colon snRNA + snATAC).
#
# PURPOSE
#   Walk the project tree, classify every file, parse contrast and cell type
#   tokens out of filenames, probe the schema of each interesting table, and
#   write a report. Nothing is interpreted here and nothing is modified. The
#   point is to see the data before deciding anything about it.
#
# WRITES
#   Only into OUT_DIR. Never touches anything under ROOT.
#
# DESIGN NOTES / WHY THINGS ARE THE WAY THEY ARE
#   - Filenames in this tree are inconsistent (spaces, typos, case, CBB vs CBD).
#     Every token goes through a normalizer rather than being matched literally.
#   - `gene-peak-matrix-peaks.txt` appears ~45 times with the same basename and
#     DIFFERENT contents depending on directory. Files are therefore always
#     identified by full path, never basename, and md5 is used to find the
#     genuine duplicates.
#   - Peak coordinates are parsed with rsplit so that scaffold names containing
#     underscores or periods, and negative TSS distances, both survive.
#   - Warn-and-continue throughout. A malformed file is reported, not fatal.
#
# Jake Lehle, 2026-10
# =============================================================================

# %%
# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

ROOT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome"
OUT_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/00_inventory"

# How many data rows to pull when probing a table's schema.
PROBE_N_ROWS = 2000

# Probe every table, or only one representative per (contrast, cell type, kind)?
# Full probing of ~200 tables takes a few minutes. Start with True.
PROBE_ALL_TABLES = True

# md5 every file to find duplicates. The FASTQs are large, so they are excluded
# from hashing by default and only inventoried by name and size.
HASH_FILES = True
HASH_SKIP_EXTENSIONS = {".gz", ".loom", ".pdf", ".jpeg", ".jpg", ".png"}
HASH_MAX_BYTES = 200 * 1024 * 1024

# Expected ATAC read layout for cellranger-arc. Used only to flag deviations.
ARC_ATAC_READS = {"R1", "R2", "R3", "I1"}
ARC_GEX_READS = {"R1", "R2", "I1", "I2"}

# Bonferroni back-calculation: how many rows to use when estimating the size of
# the tested peak universe from p_val_adj / p_val.
PEAK_UNIVERSE_PROBE_ROWS = 50

VERBOSE = True


# %%
# -----------------------------------------------------------------------------
# Cell 2: imports and small utilities
# -----------------------------------------------------------------------------

import os
import re
import csv
import sys
import json
import gzip
import hashlib
import warnings
from collections import defaultdict, Counter

import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")
os.makedirs(OUT_DIR, exist_ok=True)

WARNINGS = []


def warn(msg):
    """Record a warning and keep going. Never raise."""
    WARNINGS.append(msg)
    if VERBOSE:
        print("  [warn] " + msg)


def say(msg=""):
    if VERBOSE:
        print(msg)


def md5_of(path, max_bytes=HASH_MAX_BYTES):
    """md5 a file, or return None if skipped or unreadable."""
    ext = os.path.splitext(path)[1].lower()
    if ext in HASH_SKIP_EXTENSIONS:
        return None
    try:
        if os.path.getsize(path) > max_bytes:
            return None
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception as exc:
        warn("md5 failed on %s (%s)" % (path, exc))
        return None


def human(nbytes):
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if nbytes < 1024:
            return "%.1f%s" % (nbytes, unit)
        nbytes /= 1024.0
    return "%.1fPB" % nbytes


# %%
# -----------------------------------------------------------------------------
# Cell 3: token normalizers
#
# The tree has real inconsistencies that will silently break naive matching:
#
#   "Fibroblast S2a"      vs  "FibroblastS2a"        (space)
#   "GlutamtergicNeurons" vs  "GlutamatergicNeurons" (missing 'a' in filenames)
#   "DA_VEHvsTHCCBB_*"    vs  "DA_VEHvsTHCCBD_*"     (CBB typo in one directory)
#   "ATAC_results"        vs  "ATAC-Results"         (colon vs brain)
#   "Bcells"              vs  "B-cells"              (filename vs cluster_id)
#
# These maps are deliberately explicit rather than fuzzy. A token that does not
# normalize is reported as unknown rather than guessed at, so new spellings
# surface in the report instead of being silently absorbed.
# -----------------------------------------------------------------------------

CONTRAST_ALIASES = {
    # colon
    "thccbdvspre": "THCCBD_vs_Pre",
    "thccbdvspreinf": "THCCBD_vs_Pre",
    "vehvspre": "VEH_vs_Pre",
    "vehvspreinf": "VEH_vs_Pre",
    "vehvsthccbd": "VEH_vs_THCCBD",
    "vehvsthccbb": "VEH_vs_THCCBD",   # known typo, same comparison
    # brain
    "thccbdvsveh": "THCCBD_vs_VEH",
}

CELLTYPE_ALIASES = {
    # --- colon ---
    "bcells": "B_cells",
    "bcell": "B_cells",
    "endothelialcells1": "Endothelial_1",
    "endothelialcells2": "Endothelial_2",
    "entericgliacells": "Enteric_glia",
    "enterochromaffincells": "Enterochromaffin",
    "enterocyte1": "Enterocyte_1",
    "enterocyte2": "Enterocyte_2",
    "enterocyte3": "Enterocyte_3",
    "entercyte1": "Enterocyte_1",     # seen in an xlsx filename
    "fibroblasts1a": "Fibroblast_S1a",
    "fibroblasts1b": "Fibroblast_S1b",
    "fibroblasts1c": "Fibroblast_S1c",
    "fibroblasts2a": "Fibroblast_S2a",
    "fibroblasts2b": "Fibroblast_S2b",
    "macrophage": "Macrophage",
    "macrophages": "Macrophage",
    "nkcells": "NK_cells",
    "plasmacells": "Plasma_cells",
    "smoothmusclecells": "Smooth_muscle",
    "tcells": "T_cells",
    "unknown": "Unknown",
    # --- brain ---
    "astrocytes": "Astrocytes",
    "d1mm": "D1M_M",
    "d1ms": "D1M_S",
    "d2mm": "D2M_M",
    "d2ms": "D2M_S",
    "gabaergicneurons1": "GABAergic_1",
    "gabaergicneurons2": "GABAergic_2",
    "glutamtergicneurons": "Glutamatergic",     # typo as it appears in files
    "glutamatergicneurons": "Glutamatergic",    # correct spelling, dir names
    "interneuron": "Interneuron",
    "microglia1": "Microglia_1",
    "microglia2": "Microglia_2",
    "microglia3": "Microglia_3",
    "oligodendrocytes": "Oligodendrocytes",
    "opc": "OPC",
}


def _squash(token):
    """Lowercase and drop everything that is not a letter or digit."""
    return re.sub(r"[^a-z0-9]", "", str(token).lower())


def norm_contrast(token):
    key = _squash(token)
    return CONTRAST_ALIASES.get(key, "UNKNOWN:" + str(token))


def norm_celltype(token):
    key = _squash(token)
    return CELLTYPE_ALIASES.get(key, "UNKNOWN:" + str(token))


# %%
# -----------------------------------------------------------------------------
# Cell 4: filename classification
#
# File kinds we expect, and what each one actually holds:
#
#   da_peaks       DA_<contrast>_<celltype>.txt in ATAC_DifAccePerGroup*
#                  Signac FindMarkers() peak-level differential accessibility.
#                  Row names are peaks. This is the only file with real
#                  accessibility statistics.
#
#   assigned_peak  *_ASSIGNED_all-peak.txt
#                  One peak per line followed by one or more gene symbols.
#
#   assigned_gene  *_ASSIGNED_all-gene.txt
#                  One gene per line followed by one or more peaks.
#
#   gpm            gene-peak-matrix-peaks.txt
#                  Two different shapes share this basename:
#                    - in ALLGENES: wide, genes x cell types, cells are peak
#                      lists or NP
#                    - in a per-cell-type dir: long, gene-level muscat/edgeR
#                      RNA statistics with a joined peak column
#                  Shape is detected in Cell 7, not assumed from the path.
# -----------------------------------------------------------------------------

FASTQ_RE = re.compile(
    r"^(?P<sample>.+)_S(?P<snum>\d+)_L(?P<lane>\d{3})_"
    r"(?P<read>[RI][123])_(?P<chunk>\d{3})\.fastq\.gz$"
)

DA_RE = re.compile(r"^DA_(?P<contrast>[^_]+)_(?P<celltype>.+)\.txt$", re.IGNORECASE)

ASSIGNED_RE = re.compile(
    r"^DA_(?P<contrast>[^_]+)_(?P<celltype>.+)_ASSIGNED_all-(?P<orient>peak|gene)\.txt$",
    re.IGNORECASE,
)


def classify(path, root):
    """Return a dict describing one file. Never raises."""
    rel = os.path.relpath(path, root)
    parts = rel.split(os.sep)
    base = parts[-1]
    lower = base.lower()

    rec = {
        "rel_path": rel,
        "abs_path": path,
        "basename": base,
        "depth": len(parts) - 1,
        "parent_dir": parts[-2] if len(parts) > 1 else "",
        "top_dir": parts[0] if len(parts) > 1 else "<root>",
        "ext": os.path.splitext(base)[1].lower(),
        "size_bytes": None,
        "kind": "other",
        "tissue": None,
        "contrast_raw": None,
        "contrast": None,
        "celltype_raw": None,
        "celltype": None,
        "in_allgenes": "ALLGENES" in parts,
        "md5": None,
    }

    try:
        rec["size_bytes"] = os.path.getsize(path)
    except Exception as exc:
        warn("stat failed on %s (%s)" % (rel, exc))

    # tissue from the top-level directory name, case-insensitively
    top_l = rec["top_dir"].lower()
    if "brain" in top_l:
        rec["tissue"] = "Brain"
    elif "colon" in top_l:
        rec["tissue"] = "Colon"
    elif "raw data" in top_l:
        rec["tissue"] = "Brain" if "brain" in rel.lower() else (
            "Colon" if "colon" in rel.lower() else None
        )

    # ---- FASTQ ----
    m = FASTQ_RE.match(base)
    if m:
        rec["kind"] = "fastq"
        rec.update({
            "fq_sample": m.group("sample"),
            "fq_snum": int(m.group("snum")),
            "fq_lane": "L" + m.group("lane"),
            "fq_read": m.group("read"),
            "fq_chunk": m.group("chunk"),
        })
        s = m.group("sample")
        # arcATAC / arcGEX for brain; _NA/_PA (ATAC) and _NG/_PG (GEX) for colon
        if re.search(r"arcATAC$", s):
            rec["fq_modality"] = "ATAC"
        elif re.search(r"arcGEX$", s):
            rec["fq_modality"] = "GEX"
        elif re.search(r"_[NP]A$", s):
            rec["fq_modality"] = "ATAC"
        elif re.search(r"_[NP]G$", s):
            rec["fq_modality"] = "GEX"
        else:
            rec["fq_modality"] = "UNKNOWN"
            warn("could not infer modality for FASTQ sample %s" % s)
        anim = re.search(r"(40\d{3})", s)
        rec["fq_animal"] = anim.group(1) if anim else None
        # colon libraries carry an N / P timepoint letter; brain does not
        tp = re.search(r"_([NP])[AG]$", s)
        rec["fq_timepoint"] = {"N": "Necropsy", "P": "Pre"}.get(
            tp.group(1)) if tp else None
        return rec

    # ---- ASSIGNED mapping files (check before the generic DA_ pattern) ----
    m = ASSIGNED_RE.match(base)
    if m:
        rec["kind"] = "assigned_" + m.group("orient").lower()
        rec["contrast_raw"] = m.group("contrast")
        rec["celltype_raw"] = m.group("celltype")
        rec["contrast"] = norm_contrast(m.group("contrast"))
        rec["celltype"] = norm_celltype(m.group("celltype"))
        return rec

    # ---- gene-peak matrix (shape determined later) ----
    if lower == "gene-peak-matrix-peaks.txt":
        rec["kind"] = "gpm"
        # the cell type is the containing directory, except in ALLGENES
        if not rec["in_allgenes"] and rec["parent_dir"]:
            rec["celltype_raw"] = rec["parent_dir"]
            rec["celltype"] = norm_celltype(rec["parent_dir"])
        # contrast is whichever path component normalizes to one
        for p in parts:
            c = norm_contrast(p)
            if not c.startswith("UNKNOWN"):
                rec["contrast_raw"] = p
                rec["contrast"] = c
                break
        return rec

    # ---- FindMarkers peak-level DA ----
    m = DA_RE.match(base)
    if m:
        rec["kind"] = "da_peaks"
        rec["contrast_raw"] = m.group("contrast")
        rec["celltype_raw"] = m.group("celltype")
        rec["contrast"] = norm_contrast(m.group("contrast"))
        rec["celltype"] = norm_celltype(m.group("celltype"))
        return rec

    # ---- everything else ----
    if rec["ext"] in (".xlsx", ".xls"):
        rec["kind"] = "excel"
    elif rec["ext"] == ".loom":
        rec["kind"] = "loom"
    elif rec["ext"] == ".pdf":
        rec["kind"] = "pdf"
    elif rec["ext"] in (".jpeg", ".jpg", ".png", ".tif", ".tiff"):
        rec["kind"] = "image"
    elif rec["ext"] == ".txt":
        rec["kind"] = "text_note"

    return rec


# %%
# -----------------------------------------------------------------------------
# Cell 5: walk the tree
# -----------------------------------------------------------------------------

say("=" * 78)
say("MULTIOME_00_inventory")
say("ROOT: %s" % ROOT)
say("=" * 78)

if not os.path.isdir(ROOT):
    sys.exit("ROOT does not exist or is not a directory: %s" % ROOT)

records = []
n_dirs = 0
for dirpath, dirnames, filenames in os.walk(ROOT):
    n_dirs += 1
    dirnames.sort()
    for fn in sorted(filenames):
        records.append(classify(os.path.join(dirpath, fn), ROOT))

files = pd.DataFrame(records)
say("\nWalked %d directories, %d files, %s total"
    % (n_dirs, len(files), human(files["size_bytes"].fillna(0).sum())))

if HASH_FILES:
    say("\nHashing non-FASTQ files to find duplicates ...")
    mask = files["kind"] != "fastq"
    files.loc[mask, "md5"] = files.loc[mask, "abs_path"].map(md5_of)

say("\nFile kinds:")
for k, n in files["kind"].value_counts().items():
    say("  %-16s %5d" % (k, n))


# %%
# -----------------------------------------------------------------------------
# Cell 6: coverage matrix
#
# Which (tissue, contrast, cell type) combinations have which file kinds. This
# is the single most useful table in the report: it shows up front where the
# tree is incomplete, so missing files become a known fact rather than an empty
# panel discovered three scripts later.
# -----------------------------------------------------------------------------

tabular = files[files["kind"].isin(
    ["da_peaks", "assigned_peak", "assigned_gene", "gpm"])].copy()

# Treat the ALLGENES copies as the canonical ASSIGNED set. The per-cell-type
# subdirectories hold same-named duplicates and are an incomplete subset.
tabular["scope"] = np.where(tabular["in_allgenes"], "ALLGENES", "per_celltype")

cov = (tabular[tabular["celltype"].notna()]
       .pivot_table(index=["tissue", "contrast", "celltype"],
                    columns="kind",
                    values="rel_path",
                    aggfunc="count",
                    fill_value=0)
       .reset_index())

say("\nCoverage by tissue / contrast / cell type:")
say(cov.to_string(index=False))

# unknown tokens are reported, never guessed
unknown_ct = sorted(set(
    t for t in tabular["celltype"].dropna().unique() if str(t).startswith("UNKNOWN")))
unknown_co = sorted(set(
    t for t in tabular["contrast"].dropna().unique() if str(t).startswith("UNKNOWN")))
for t in unknown_ct:
    warn("unrecognized cell type token: %s  (add to CELLTYPE_ALIASES)" % t)
for t in unknown_co:
    warn("unrecognized contrast token: %s  (add to CONTRAST_ALIASES)" % t)

# cell type set per contrast, and the symmetric difference between contrasts
ct_sets = {}
for (tis, con), grp in tabular.groupby(["tissue", "contrast"]):
    ct_sets[(tis, con)] = set(grp["celltype"].dropna().unique())

say("\nCell type set size per contrast:")
for k in sorted(ct_sets):
    say("  %-8s %-18s %3d" % (k[0], k[1], len(ct_sets[k])))

say("\nCell types NOT shared across all contrasts within a tissue:")
for tis in sorted(set(k[0] for k in ct_sets if k[0])):
    keys = [k for k in ct_sets if k[0] == tis]
    if len(keys) < 2:
        continue
    common = set.intersection(*[ct_sets[k] for k in keys])
    union = set.union(*[ct_sets[k] for k in keys])
    for ct in sorted(union - common):
        present = [k[1] for k in keys if ct in ct_sets[k]]
        say("  %-8s %-18s present only in: %s" % (tis, ct, ", ".join(sorted(present))))


# %%
# -----------------------------------------------------------------------------
# Cell 7: probe table schemas
#
# Sniffs delimiter, header presence, column names, dtypes and row counts. For
# the gene-peak matrix it also decides whether the file is the wide cell-type
# matrix or the long per-cell-type statistics table, since both share a
# basename.
# -----------------------------------------------------------------------------

def sniff_delimiter(path, n_lines=20):
    """Count candidate delimiters across the first few lines."""
    try:
        with open(path, "r", errors="replace") as fh:
            head = [fh.readline() for _ in range(n_lines)]
    except Exception as exc:
        warn("could not read %s (%s)" % (path, exc))
        return None
    head = [h for h in head if h.strip()]
    if not head:
        return None
    scores = {}
    for d, name in [("\t", "tab"), (",", "comma")]:
        counts = [h.count(d) for h in head]
        # a real delimiter appears consistently and more than zero times
        if min(counts) > 0 and len(set(counts)) <= 3:
            scores[name] = np.mean(counts)
    if scores:
        return max(scores, key=scores.get)
    # fall back to whitespace if lines split into a stable number of fields
    nf = [len(h.split()) for h in head]
    if min(nf) >= 2:
        return "whitespace"
    return None


DELIM_MAP = {"tab": "\t", "comma": ",", "whitespace": None}


def probe_table(path, kind):
    """Return a schema dict for one tabular file. Never raises."""
    out = {
        "rel_path": os.path.relpath(path, ROOT),
        "kind": kind,
        "delimiter": None,
        "n_cols": None,
        "columns": None,
        "has_header": None,
        "n_rows_probed": None,
        "n_rows_total": None,
        "shape_call": None,
        "notes": [],
    }

    delim = sniff_delimiter(path)
    out["delimiter"] = delim
    if delim is None:
        warn("no delimiter detected in %s" % out["rel_path"])
        return out

    sep = DELIM_MAP[delim]
    read_kw = dict(nrows=PROBE_N_ROWS, dtype=str, engine="python")
    if sep is None:
        read_kw["sep"] = r"\s+"
    else:
        read_kw["sep"] = sep

    # total line count, cheaply
    try:
        with open(path, "rb") as fh:
            out["n_rows_total"] = sum(1 for _ in fh)
    except Exception as exc:
        warn("line count failed on %s (%s)" % (out["rel_path"], exc))

    # ASSIGNED files are ragged (variable fields per line) so pandas is wrong
    # for them. Handle separately.
    if kind in ("assigned_peak", "assigned_gene"):
        try:
            widths, extras = [], []
            with open(path, "r", errors="replace") as fh:
                for i, line in enumerate(fh):
                    if i >= PROBE_N_ROWS:
                        break
                    f = line.split()
                    if not f:
                        continue
                    widths.append(len(f))
                    extras.append(len(f) - 1)
            out["n_rows_probed"] = len(widths)
            out["n_cols"] = "ragged %d-%d" % (min(widths), max(widths)) if widths else None
            out["shape_call"] = "ragged_mapping"
            out["notes"].append(
                "fields per line: min %d, max %d, median %.1f"
                % (min(widths), max(widths), float(np.median(widths))) if widths else "empty")
            out["notes"].append(
                "multi-mapping lines (>1 partner): %d of %d"
                % (sum(1 for e in extras if e > 1), len(extras)))
            out["columns"] = ["key", "partners..."]
        except Exception as exc:
            warn("ragged probe failed on %s (%s)" % (out["rel_path"], exc))
        return out

    try:
        df = pd.read_csv(path, **read_kw)
    except Exception as exc:
        warn("pandas failed on %s (%s)" % (out["rel_path"], exc))
        return out

    out["n_cols"] = df.shape[1]
    out["n_rows_probed"] = df.shape[0]
    out["columns"] = list(df.columns)

    # Header detection: FindMarkers output written by write.table with
    # row.names=TRUE and col.names=TRUE has one fewer header field than data
    # field, which pandas resolves by making the peak column the index.
    first = str(df.columns[0])
    out["has_header"] = not bool(re.match(r"^chr[\w.]+[-:]\d+[-]\d+$", first))

    # ---- decide what a gene-peak matrix actually is ----
    if kind == "gpm":
        cols_l = [str(c).lower() for c in df.columns]
        stat_cols = {"logfc", "logcpm", "p_val", "p_adj.loc", "p_adj.glb",
                     "cluster_id", "contrast"}
        n_da_cols = sum(1 for c in df.columns if str(c).startswith("DA_"))
        if stat_cols.intersection(cols_l):
            out["shape_call"] = "long_gene_stats"
            out["notes"].append(
                "gene-level RNA statistics (muscat/edgeR signature), "
                "peak column is a joined annotation")
        elif n_da_cols >= 2:
            out["shape_call"] = "wide_celltype_matrix"
            out["notes"].append("%d DA_* cell type columns" % n_da_cols)
        elif n_da_cols == 1:
            out["shape_call"] = "long_single_celltype"
        else:
            out["shape_call"] = "unrecognized"
            warn("could not classify gpm shape: %s" % out["rel_path"])

    if kind == "da_peaks":
        cols_l = [str(c).lower() for c in df.columns]
        if "avg_log2fc" in cols_l and "pct.1" in cols_l:
            out["shape_call"] = "findmarkers"
        elif "logfc" in cols_l and "p_adj.loc" in cols_l:
            out["shape_call"] = "muscat_pbds"
            out["notes"].append(
                "muscat columns in a DA_ file: confirm whether this is peak "
                "level or gene level before using it")
        else:
            out["shape_call"] = "unrecognized"
            warn("unrecognized DA_ schema: %s -> %s"
                 % (out["rel_path"], list(df.columns)))

    # ---- NP token and peak column presence ----
    for c in df.columns:
        vals = df[c].dropna().astype(str)
        if len(vals) and (vals == "NP").mean() > 0.01:
            out["notes"].append("column '%s' is %.1f%% NP"
                                % (c, 100 * (vals == "NP").mean()))

    return out


probe_targets = tabular.copy()
if not PROBE_ALL_TABLES:
    probe_targets = probe_targets.drop_duplicates(
        subset=["tissue", "contrast", "kind", "in_allgenes"])

say("\nProbing %d tables ..." % len(probe_targets))
schemas = [probe_table(r.abs_path, r.kind) for r in probe_targets.itertuples()]
schema_df = pd.DataFrame(schemas)

say("\nDistinct schemas found:")
sig = schema_df.copy()
sig["colsig"] = sig["columns"].map(lambda c: " | ".join(map(str, c)) if isinstance(c, list) else str(c))
for (kind, shape, colsig), grp in sig.groupby(["kind", "shape_call", "colsig"], dropna=False):
    say("\n  kind=%s  shape=%s  n_files=%d" % (kind, shape, len(grp)))
    say("    columns: %s" % colsig)
    say("    example: %s" % grp["rel_path"].iloc[0])
    say("    rows (total incl. header): %s"
        % sorted(set(grp["n_rows_total"].dropna().astype(int)))[:8])


# %%
# -----------------------------------------------------------------------------
# Cell 8: duplicate detection
#
# The per-cell-type subdirectories appear to hold byte-identical copies of the
# ALLGENES files. Confirming that means the pipeline can read ALLGENES only and
# ignore roughly half the tree.
# -----------------------------------------------------------------------------

if HASH_FILES:
    dup = (files[files["md5"].notna()]
           .groupby("md5")
           .agg(n=("rel_path", "size"), paths=("rel_path", list))
           .query("n > 1")
           .sort_values("n", ascending=False))
    say("\n%d md5 groups with more than one copy" % len(dup))
    for h, row in dup.head(25).iterrows():
        say("  %dx  %s" % (row["n"], row["paths"][0]))
        for p in row["paths"][1:4]:
            say("       = %s" % p)
        if row["n"] > 4:
            say("       ... and %d more" % (row["n"] - 4))

    # specifically: is every per-cell-type ASSIGNED file a copy of an ALLGENES one?
    asg = files[files["kind"].isin(["assigned_peak", "assigned_gene"])]
    ag_hashes = set(asg[asg["in_allgenes"]]["md5"].dropna())
    pc = asg[~asg["in_allgenes"]]
    matched = pc["md5"].isin(ag_hashes).sum()
    say("\nper-cell-type ASSIGNED files identical to an ALLGENES copy: %d of %d"
        % (matched, len(pc)))
    if len(pc) and matched < len(pc):
        for p in pc[~pc["md5"].isin(ag_hashes)]["rel_path"].head(10):
            warn("per-cell-type ASSIGNED file is NOT a copy of ALLGENES: %s" % p)


# %%
# -----------------------------------------------------------------------------
# Cell 9: peak coordinate parsing and chromosome inventory
#
# Validates the coordinate parser against real data and inventories contig
# names. Scaffold names can contain underscores and periods, and TSS distances
# can be negative, so rsplit is used rather than split.
# -----------------------------------------------------------------------------

PEAK_WITH_DIST_RE = re.compile(r"^(?P<peak>.+?)_(?P<dist>-?\d+)$")


def parse_peak(token):
    """
    'chr1-56999-57879_-25'  -> ('chr1', 56999, 57879, -25)
    'chr1-56999-57879'      -> ('chr1', 56999, 57879, None)
    Returns None if the token does not parse.
    """
    token = str(token).strip()
    if not token or token == "NP":
        return None
    dist = None
    m = PEAK_WITH_DIST_RE.match(token)
    if m:
        token = m.group("peak")
        dist = int(m.group("dist"))
    # chrom may itself contain '-'; take the last two fields as start/end
    bits = token.rsplit("-", 2)
    if len(bits) != 3:
        return None
    chrom, start, end = bits
    try:
        return (chrom, int(start), int(end), dist)
    except ValueError:
        return None


def inventory_peaks(path, kind, max_rows=PROBE_N_ROWS):
    """Pull peak tokens out of one file and summarize them."""
    chroms, widths, dists, n_bad, n_tok = Counter(), [], [], 0, 0
    try:
        with open(path, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= max_rows:
                    break
                for tok in re.split(r"[\t;]| {2,}|\s", line.strip()):
                    if not tok.startswith("chr") and not tok.startswith("NW_") \
                       and not tok.startswith("NC_"):
                        continue
                    n_tok += 1
                    p = parse_peak(tok)
                    if p is None:
                        n_bad += 1
                        continue
                    chroms[p[0]] += 1
                    widths.append(p[2] - p[1])
                    if p[3] is not None:
                        dists.append(p[3])
    except Exception as exc:
        warn("peak inventory failed on %s (%s)" % (path, exc))
    return {
        "rel_path": os.path.relpath(path, ROOT),
        "kind": kind,
        "n_tokens": n_tok,
        "n_unparsed": n_bad,
        "n_contigs": len(chroms),
        "contigs": sorted(chroms),
        "width_min": min(widths) if widths else None,
        "width_median": float(np.median(widths)) if widths else None,
        "width_max": max(widths) if widths else None,
        "dist_min": min(dists) if dists else None,
        "dist_max": max(dists) if dists else None,
        "dist_n_negative": sum(1 for d in dists if d < 0),
        "dist_n_total": len(dists),
    }


peak_targets = tabular.drop_duplicates(subset=["tissue", "contrast", "kind", "in_allgenes"])
say("\nInventorying peak coordinates across %d representative files ..." % len(peak_targets))
peak_inv = pd.DataFrame([inventory_peaks(r.abs_path, r.kind)
                         for r in peak_targets.itertuples()])

all_contigs = sorted(set(c for lst in peak_inv["contigs"] for c in lst))
say("\n%d distinct contig names seen: %s" % (len(all_contigs), all_contigs))
non_chr = [c for c in all_contigs if not re.match(r"^chr([0-9]+|[XYM]|MT)$", c)]
if non_chr:
    warn("non-standard contig names present: %s" % non_chr)

say("\nPeak width and TSS distance ranges:")
say(peak_inv[["rel_path", "kind", "n_tokens", "n_unparsed", "width_min",
              "width_median", "width_max", "dist_min", "dist_max",
              "dist_n_negative", "dist_n_total"]].to_string(index=False))

bad = peak_inv[peak_inv["n_unparsed"] > 0]
for r in bad.itertuples():
    warn("%d peak tokens failed to parse in %s" % (r.n_unparsed, r.rel_path))


# %%
# -----------------------------------------------------------------------------
# Cell 10: recover the tested peak universe
#
# Signac FindMarkers() adjusts with Bonferroni, so p_val_adj / p_val recovers
# the number of tests, which is the size of the peak set that went into the
# comparison. That gives a denominator for cross-cell-type comparison without
# needing the peak BED.
#
# Also flags the pseudocount-floor rows: pct.1 or pct.2 exactly zero produces an
# avg_log2FC set by the pseudocount rather than by effect size, and those rows
# will dominate any volcano or ranked bar plot if left in.
# -----------------------------------------------------------------------------

def peak_universe(path):
    out = {"rel_path": os.path.relpath(path, ROOT)}
    try:
        df = pd.read_csv(path, sep="\t", index_col=0, nrows=PEAK_UNIVERSE_PROBE_ROWS)
    except Exception as exc:
        warn("universe probe failed on %s (%s)" % (out["rel_path"], exc))
        return out
    cols = {str(c).lower(): c for c in df.columns}
    if "p_val" not in cols or "p_val_adj" not in cols:
        return out
    p, padj = df[cols["p_val"]], df[cols["p_val_adj"]]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (padj / p).replace([np.inf, -np.inf], np.nan).dropna()
    # Bonferroni caps at 1.0, so only use rows well below the cap
    ratio = ratio[padj.reindex(ratio.index) < 0.9]
    if len(ratio):
        out["n_tests_median"] = float(np.median(ratio))
        out["n_tests_min"] = float(ratio.min())
        out["n_tests_max"] = float(ratio.max())
        out["n_tests_rows_used"] = int(len(ratio))

    # pseudocount floor
    if "pct.1" in cols and "pct.2" in cols:
        z = (df[cols["pct.1"]] == 0) | (df[cols["pct.2"]] == 0)
        out["n_zero_pct_rows"] = int(z.sum())
        if "avg_log2fc" in cols and z.any():
            out["floor_log2fc_min"] = float(df.loc[z, cols["avg_log2fc"]].min())
            out["floor_log2fc_max"] = float(df.loc[z, cols["avg_log2fc"]].max())
            nz = ~z
            if nz.any():
                out["real_log2fc_min"] = float(df.loc[nz, cols["avg_log2fc"]].min())
                out["real_log2fc_max"] = float(df.loc[nz, cols["avg_log2fc"]].max())
    return out


da_files = files[files["kind"] == "da_peaks"]
say("\nRecovering peak universe from %d FindMarkers files ..." % len(da_files))
universe = pd.DataFrame([peak_universe(p) for p in da_files["abs_path"]])
universe = universe.merge(
    da_files[["rel_path", "tissue", "contrast", "celltype"]], on="rel_path", how="left")

cols_show = [c for c in ["tissue", "contrast", "celltype", "n_tests_median",
                         "n_zero_pct_rows", "floor_log2fc_min", "real_log2fc_min",
                         "real_log2fc_max"] if c in universe.columns]
say("\nPeak universe and pseudocount-floor summary:")
say(universe[cols_show].to_string(index=False))

if "n_tests_median" in universe.columns:
    vals = universe["n_tests_median"].dropna()
    if len(vals):
        say("\nEstimated tested peaks: median %.0f, range %.0f to %.0f across files"
            % (np.median(vals), vals.min(), vals.max()))
        if vals.max() - vals.min() > 1000:
            warn("peak universe size varies by more than 1000 across cell types; "
                 "peaks were probably filtered per cell type, so counts are not "
                 "directly comparable without the per-cell-type denominator")


# %%
# -----------------------------------------------------------------------------
# Cell 11: FASTQ inventory and cellranger-arc libraries.csv drafts
#
# One cellranger-arc run per (tissue, animal, timepoint). Writes a draft
# libraries.csv per run into OUT_DIR for review. Nothing is launched here.
# -----------------------------------------------------------------------------

fq = files[files["kind"] == "fastq"].copy()
say("\nFASTQ files: %d, %s total" % (len(fq), human(fq["size_bytes"].fillna(0).sum())))

if len(fq):
    fq["run_key"] = (fq["tissue"].fillna("?") + "_" + fq["fq_animal"].fillna("?")
                     + "_" + fq["fq_timepoint"].fillna("single"))

    lib = (fq.groupby(["run_key", "tissue", "fq_animal", "fq_timepoint",
                       "fq_modality", "fq_sample"], dropna=False)
             .agg(n_files=("rel_path", "size"),
                  lanes=("fq_lane", lambda s: ",".join(sorted(set(s)))),
                  reads=("fq_read", lambda s: ",".join(sorted(set(s)))),
                  size=("size_bytes", "sum"),
                  fastq_dir=("abs_path", lambda s: os.path.dirname(s.iloc[0])))
             .reset_index())
    say("\nLibraries per run:")
    say(lib[["run_key", "fq_modality", "fq_sample", "lanes", "reads", "n_files"]]
        .to_string(index=False))

    # read layout check
    for r in lib.itertuples():
        got = set(r.reads.split(","))
        want = ARC_ATAC_READS if r.fq_modality == "ATAC" else ARC_GEX_READS
        if got != want:
            warn("%s (%s): read files are %s, cellranger-arc expects %s"
                 % (r.fq_sample, r.fq_modality, sorted(got), sorted(want)))

    # runs missing a modality
    for rk, grp in lib.groupby("run_key"):
        mods = set(grp["fq_modality"])
        if mods != {"ATAC", "GEX"}:
            warn("run %s has modalities %s; cellranger-arc needs both"
                 % (rk, sorted(mods)))

    # draft libraries.csv per run
    lib_dir = os.path.join(OUT_DIR, "libraries_draft")
    os.makedirs(lib_dir, exist_ok=True)
    for rk, grp in lib.groupby("run_key"):
        path = os.path.join(lib_dir, "%s_libraries.csv" % rk)
        try:
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(["fastqs", "sample", "library_type"])
                for r in grp.itertuples():
                    lt = ("Chromatin Accessibility" if r.fq_modality == "ATAC"
                          else "Gene Expression")
                    w.writerow([r.fastq_dir, r.fq_sample, lt])
        except Exception as exc:
            warn("could not write %s (%s)" % (path, exc))
    say("\nWrote %d draft libraries.csv files to %s"
        % (lib["run_key"].nunique(), lib_dir))

    say("\nAnimal by timepoint coverage:")
    say(pd.crosstab([fq["tissue"], fq["fq_animal"]],
                    [fq["fq_timepoint"].fillna("single"), fq["fq_modality"]]).to_string())
else:
    lib = pd.DataFrame()
    warn("no FASTQ files matched the Illumina naming pattern")


# %%
# -----------------------------------------------------------------------------
# Cell 12: metadata, notes and loom headers
#
# Reads the small text and Excel files that are most likely to answer the open
# methods questions: which reference build, what thresholds, which group was
# ident.1. Reads loom headers only, not the matrices.
# -----------------------------------------------------------------------------

say("\n" + "-" * 78)
say("Text notes")
say("-" * 78)
for r in files[files["kind"] == "text_note"].itertuples():
    say("\n### %s" % r.rel_path)
    try:
        with open(r.abs_path, "r", errors="replace") as fh:
            txt = fh.read(8000)
        say(txt if txt.strip() else "  (empty)")
    except Exception as exc:
        warn("could not read %s (%s)" % (r.rel_path, exc))

say("\n" + "-" * 78)
say("Excel files: sheet names and headers")
say("-" * 78)
excel_summary = []
for r in files[files["kind"] == "excel"].itertuples():
    say("\n### %s" % r.rel_path)
    try:
        xl = pd.ExcelFile(r.abs_path)
        for sheet in xl.sheet_names:
            try:
                d = xl.parse(sheet, nrows=5)
                say("  [%s]  %d cols: %s" % (sheet, d.shape[1], list(d.columns)[:12]))
                excel_summary.append({
                    "rel_path": r.rel_path, "sheet": sheet,
                    "n_cols": d.shape[1], "columns": list(d.columns),
                })
            except Exception as exc:
                warn("sheet %s of %s failed (%s)" % (sheet, r.rel_path, exc))
    except Exception as exc:
        warn("could not open %s (%s)" % (r.rel_path, exc))

say("\n" + "-" * 78)
say("Loom files: header only")
say("-" * 78)
for r in files[files["kind"] == "loom"].itertuples():
    say("\n### %s  (%s)" % (r.rel_path, human(r.size_bytes or 0)))
    try:
        import h5py
        with h5py.File(r.abs_path, "r") as h5:
            if "matrix" in h5:
                say("  matrix shape: %s" % (h5["matrix"].shape,))
            for grp in ("col_attrs", "row_attrs"):
                if grp in h5:
                    say("  %s: %s" % (grp, sorted(h5[grp].keys())))
            if "layers" in h5:
                say("  layers: %s" % sorted(h5["layers"].keys()))
    except ImportError:
        warn("h5py not available in this env; skipping loom header read")
    except Exception as exc:
        warn("could not read loom header %s (%s)" % (r.rel_path, exc))

say("\nPDF reports present (read these by hand for the reference build):")
for r in files[files["kind"] == "pdf"].itertuples():
    say("  %s  (%s)" % (r.rel_path, human(r.size_bytes or 0)))

say("\nImages present:")
for r in files[files["kind"] == "image"].itertuples():
    say("  %s" % r.rel_path)


# %%
# -----------------------------------------------------------------------------
# Cell 13: write the report
# -----------------------------------------------------------------------------

def dump(df, name):
    if df is None or not len(df):
        return
    path = os.path.join(OUT_DIR, name)
    try:
        out = df.copy()
        for c in out.columns:
            if out[c].map(lambda v: isinstance(v, (list, dict, set))).any():
                out[c] = out[c].map(
                    lambda v: ";".join(map(str, v)) if isinstance(v, (list, set))
                    else json.dumps(v) if isinstance(v, dict) else v)
        out.to_csv(path, index=False, lineterminator="\n")
        say("  wrote %s  (%d rows)" % (name, len(out)))
    except Exception as exc:
        warn("could not write %s (%s)" % (name, exc))


say("\n" + "=" * 78)
say("Writing report to %s" % OUT_DIR)
say("=" * 78)

dump(files.drop(columns=["abs_path"], errors="ignore"), "file_inventory.csv")
dump(cov, "coverage_matrix.csv")
dump(schema_df, "table_schemas.csv")
dump(peak_inv, "peak_coordinate_inventory.csv")
dump(universe, "peak_universe_and_floor.csv")
dump(lib, "fastq_libraries.csv")
dump(pd.DataFrame(excel_summary), "excel_sheet_summary.csv")
dump(pd.DataFrame({"warning": WARNINGS}), "warnings.csv")

summary = {
    "root": ROOT,
    "generated": pd.Timestamp.now().isoformat(),
    "n_files": int(len(files)),
    "n_dirs": int(n_dirs),
    "total_bytes": int(files["size_bytes"].fillna(0).sum()),
    "file_kinds": {k: int(v) for k, v in files["kind"].value_counts().items()},
    "tissues": sorted(set(files["tissue"].dropna())),
    "contrasts": sorted(set(tabular["contrast"].dropna())),
    "celltypes_by_contrast": {
        "%s/%s" % k: sorted(v) for k, v in ct_sets.items()},
    "contigs": all_contigs,
    "n_warnings": len(WARNINGS),
}
try:
    with open(os.path.join(OUT_DIR, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    say("  wrote summary.json")
except Exception as exc:
    warn("could not write summary.json (%s)" % exc)

say("\n%d warnings raised. Read warnings.csv before using any of this." % len(WARNINGS))
say("Done.")
