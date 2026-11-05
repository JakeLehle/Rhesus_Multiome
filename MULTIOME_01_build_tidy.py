# %%
# =============================================================================
# MULTIOME_01_build_tidy.py
#
# Read-only parser for the THC/CBD rhesus multiome collaboration. Turns the
# vendor's scattered tables into four tidy long CSVs plus a diagnostic report.
# No figures, no filtering, no interpretation.
#
# WHAT EACH SOURCE FILE ACTUALLY CONTAINS
#   (confirmed by the vendor note at ATAC_results_Brain/ATAC-Results/Note on results.txt
#    and by the MULTIOME_00 schema probe)
#
#   DA_<contrast>_<celltype>.txt
#       Signac FindMarkers() on the peak assay. Row names are peaks.
#       THE ONLY source of accessibility statistics. Already filtered by the
#       vendor to their DA call, so this script records the statistics and does
#       not re-threshold.
#
#   ALLGENES/*_ASSIGNED_all-peak.txt
#       Peak -> gene(s). The gene side covers all Mmul genes, so this is the
#       complete peak-to-gene map for the DA peaks.
#
#   ALLGENES/*_ASSIGNED_all-gene.txt
#       Gene -> peak(s). The reciprocal view, used here only to validate.
#
#   <celltype>/*_ASSIGNED_all-*.txt
#       Same relation but restricted to genes in that cell type's tested RNA
#       universe. NOT a duplicate of the ALLGENES copy.
#
#   ALLGENES/gene-peak-matrix-peaks.txt
#       Wide, 20596 genes x cell types. Cells hold "NP" or
#       "chr-start-end_distanceToTSS" lists. THE ONLY source of TSS distance.
#
#   <celltype>/gene-peak-matrix-peaks.txt
#       Long, that cell type's tested RNA gene universe with muscat/edgeR
#       statistics (logFC here is EXPRESSION, not accessibility) plus the peak
#       annotation column. Row counts are identical across contrasts for a
#       given cell type, which confirms it is the post-filterByExpr tested set
#       rather than a significant subset.
#
# OUTPUTS (into OUT_DIR)
#   peaks_da.csv         one row per tissue/contrast/celltype/peak
#   peak_gene_map.csv    one row per tissue/contrast/celltype/peak/gene
#   genes_rna.csv        one row per tissue/contrast/celltype/gene
#   integrated.csv       the outer join of the three, nothing dropped
#   plus per-group summaries, validation tables and warnings.csv
#
# Jake Lehle, 2026-10
# =============================================================================

# %%
# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

ROOT = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome"
OUT_DIR = "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/01_tidy"

# Restrict the run, or leave empty for everything.
#   e.g. TISSUES = ["Brain"]; CONTRASTS = ["VEH_vs_THCCBD"]
TISSUES = []
CONTRASTS = []

# --- pseudocount floor ---------------------------------------------------
# A peak detected in almost no cells of one group gets an avg_log2FC set by the
# pseudocount rather than by effect size. MULTIOME_00 showed that an exact-zero
# test is too strict: non-zero-pct rows still reached |avg_log2FC| of 11.6.
# Any row whose smaller pct is below this is flagged floor_flag = True. Nothing
# is dropped; the flag exists so figures can separate them.
MIN_PCT_FLOOR = 0.02

# A second, independent floor check on effect size alone. Rows above this are
# flagged regardless of pct, because an honest single-cell ATAC log2FC does not
# reach double digits.
MAX_PLAUSIBLE_ABS_LOG2FC = 8.0

# --- region classification -----------------------------------------------
# Negative distance is upstream (vendor note). Observed window is about
# -20000 to +2000, so nothing outside that was ever annotated.
PROMOTER_BP = 1000
PROXIMAL_BP = 10000

# --- orientation of the ATAC effect size ---------------------------------
# UNRESOLVED. The FindMarkers files carry no record of which group was ident.1,
# so the sign of avg_log2FC cannot be read as "opened" or "closed" yet.
# Leave as None until Lakmini confirms. When set to the group name that was
# ident.1, an oriented column is written and the orientation is recorded in the
# output so no downstream script can forget it.
ATAC_IDENT1_GROUP = None

# Validate the two ASSIGNED orientations against each other. Slower, worth
# doing once.
VALIDATE_RECIPROCAL = True

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

# multiome_common.py must sit next to this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) or ".")
from multiome_common import (
    norm_contrast, norm_celltype, infer_tissue,
    parse_peak, classify_region, peek_gpm_header,
    DA_RE, ASSIGNED_RE, BRAIN_LABEL_WARNINGS, OBSERVED_WINDOW,
)

warnings.filterwarnings("ignore")
os.makedirs(OUT_DIR, exist_ok=True)

WARNINGS = []


def warn(msg):
    WARNINGS.append(msg)
    if VERBOSE:
        print("  [warn] " + msg)


def say(msg=""):
    if VERBOSE:
        print(msg)


def keep(tissue, contrast):
    if TISSUES and tissue not in TISSUES:
        return False
    if CONTRASTS and contrast not in CONTRASTS:
        return False
    return True


# %%
# -----------------------------------------------------------------------------
# Cell 3: discover the source files
#
# Three buckets, each keyed on (tissue, contrast, celltype):
#   da        the FindMarkers tables
#   asg_peak  the ALLGENES peak -> gene maps
#   asg_gene  the ALLGENES gene -> peak maps
# plus one wide matrix per (tissue, contrast) and one long gpm per key.
# -----------------------------------------------------------------------------

say("=" * 78)
say("MULTIOME_01_build_tidy")
say("ROOT: %s" % ROOT)
say("=" * 78)

if not os.path.isdir(ROOT):
    sys.exit("ROOT not found: %s" % ROOT)

src_da, src_asg_peak, src_asg_gene = {}, {}, {}
src_gpm_long, src_gpm_wide = {}, {}

for dirpath, dirnames, filenames in os.walk(ROOT):
    dirnames.sort()
    rel_dir = os.path.relpath(dirpath, ROOT)
    parts = [] if rel_dir == "." else rel_dir.split(os.sep)
    in_allgenes = "ALLGENES" in parts

    for fn in sorted(filenames):
        path = os.path.join(dirpath, fn)
        pp = parts + [fn]
        tissue = infer_tissue(pp)

        if fn.lower() == "gene-peak-matrix-peaks.txt":
            shape, contrast, celltypes = peek_gpm_header(path)
            if shape is None:
                warn("could not read gpm header: %s" % os.path.relpath(path, ROOT))
                continue
            if tissue is None:
                # the root-level copy has no tissue in its path; the cell type
                # tokens identify it
                tissue = "Brain" if any(
                    c in ("Astrocytes", "D1M_M", "OPC") for c in celltypes) else "Colon"
            if shape == "wide":
                src_gpm_wide.setdefault((tissue, contrast), []).append(path)
            else:
                src_gpm_long[(tissue, contrast, celltypes[0])] = path
            continue

        m = ASSIGNED_RE.match(fn)
        if m:
            if not in_allgenes:
                continue          # the DE-restricted copy; ALLGENES is the complete map
            key = (tissue, norm_contrast(m.group("contrast")),
                   norm_celltype(m.group("celltype")))
            (src_asg_peak if m.group("orient").lower() == "peak"
             else src_asg_gene)[key] = path
            continue

        m = DA_RE.match(fn)
        if m and "DifAccePerGroup" in rel_dir:
            key = (tissue, norm_contrast(m.group("contrast")),
                   norm_celltype(m.group("celltype")))
            src_da[key] = path

say("\nDiscovered:")
say("  FindMarkers DA tables      %4d" % len(src_da))
say("  ALLGENES peak->gene maps   %4d" % len(src_asg_peak))
say("  ALLGENES gene->peak maps   %4d" % len(src_asg_gene))
say("  per-celltype RNA tables    %4d" % len(src_gpm_long))
say("  wide gene-peak matrices    %4d groups" % len(src_gpm_wide))

for (tis, con), paths in sorted(src_gpm_wide.items()):
    if len(paths) > 1:
        say("  %s / %s: %d wide copies, using the ALLGENES one"
            % (tis, con, len(paths)))

# unknown tokens
for key in list(src_da) + list(src_asg_peak):
    for part in key[1:]:
        if str(part).startswith("UNKNOWN"):
            warn("unrecognized token in key %s (extend the alias maps)" % (key,))

# brain label discrepancies, reported not applied
if any(k[0] == "Brain" for k in src_da):
    say("\nBrain cluster label discrepancies (filename token -> "
        "ATAC_all_clusters.xlsx sheet name):")
    for tok, sheet in BRAIN_LABEL_WARNINGS.items():
        say("    %-24s -> %s" % (tok, sheet))
    warn("brain filename tokens disagree with ATAC_all_clusters.xlsx sheet names; "
         "the cluster written as GlutamtergicNeurons maps to a sheet named "
         "'GABAergic CCK+ Neurons'. Verify against Brain_cluster_markers.xlsx "
         "('Initial name' vs 'New name') before any figure is labeled.")


# %%
# -----------------------------------------------------------------------------
# Cell 4: FindMarkers peak-level DA
#
# These were written by R's write.table with row names, so the header has one
# fewer field than the data rows. Handled explicitly rather than relying on
# pandas guessing.
# -----------------------------------------------------------------------------

DA_COLS = ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]


def read_da(path):
    """Return a DataFrame indexed by peak id, or None."""
    try:
        with open(path, "r", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            first = fh.readline().rstrip("\n").split("\t")
    except Exception as exc:
        warn("could not open %s (%s)" % (os.path.relpath(path, ROOT), exc))
        return None

    if not first or first == [""]:
        return pd.DataFrame(columns=DA_COLS)          # header only

    idx = 0 if len(first) == len(header) + 1 else None
    try:
        df = pd.read_csv(path, sep="\t", index_col=idx)
    except Exception as exc:
        warn("pandas failed on %s (%s)" % (os.path.relpath(path, ROOT), exc))
        return None

    if idx is None:
        # a named first column; find whichever one holds peak coordinates
        peak_col = next((c for c in df.columns
                         if df[c].astype(str).str.match(r"^chr[\w.]+-\d+-\d+$").mean() > 0.9),
                        None)
        if peak_col is None:
            warn("no peak column found in %s" % os.path.relpath(path, ROOT))
            return None
        df = df.set_index(peak_col)

    df.index.name = "peak_id"
    missing = [c for c in DA_COLS if c not in df.columns]
    if missing:
        warn("%s lacks expected columns %s" % (os.path.relpath(path, ROOT), missing))
    return df


say("\n" + "-" * 78)
say("Loading %d FindMarkers tables" % len(src_da))
say("-" * 78)

rows = []
for (tis, con, ct), path in sorted(src_da.items()):
    if not keep(tis, con):
        continue
    df = read_da(path)
    if df is None:
        continue
    if not len(df):
        warn("%s / %s / %s: DA table is empty (header only, no peaks)" % (tis, con, ct))
        continue

    for peak_id, r in df.iterrows():
        p = parse_peak(peak_id)
        if p is None:
            warn("unparsable peak '%s' in %s / %s / %s" % (peak_id, tis, con, ct))
            continue
        _, chrom, start, end, _ = p
        pct1 = r.get("pct.1", np.nan)
        pct2 = r.get("pct.2", np.nan)
        lfc = r.get("avg_log2FC", np.nan)
        min_pct = np.nanmin([pct1, pct2]) if not (pd.isna(pct1) and pd.isna(pct2)) else np.nan
        rows.append({
            "tissue": tis, "contrast": con, "celltype": ct,
            "peak_id": p[0], "chrom": chrom, "start": start, "end": end,
            "width": end - start,
            "atac_log2FC": lfc,
            "atac_p_val": r.get("p_val", np.nan),
            "atac_p_val_adj": r.get("p_val_adj", np.nan),
            "pct_1": pct1, "pct_2": pct2, "min_pct": min_pct,
            "floor_flag": bool(
                (not pd.isna(min_pct) and min_pct < MIN_PCT_FLOOR)
                or (not pd.isna(lfc) and abs(lfc) > MAX_PLAUSIBLE_ABS_LOG2FC)),
            "atac_direction_raw": ("up_in_pct1" if (not pd.isna(lfc) and lfc > 0)
                                   else "down_in_pct1" if not pd.isna(lfc) else None),
        })

peaks_da = pd.DataFrame(rows)
say("\npeaks_da: %d rows, %d unique peaks"
    % (len(peaks_da), peaks_da["peak_id"].nunique() if len(peaks_da) else 0))

if len(peaks_da):
    say("floor-flagged rows: %d of %d (%.1f%%)"
        % (peaks_da["floor_flag"].sum(), len(peaks_da),
           100 * peaks_da["floor_flag"].mean()))

    # recover the tested peak universe from the Bonferroni ratio, per group
    uni = []
    for (tis, con, ct), g in peaks_da.groupby(["tissue", "contrast", "celltype"]):
        r = (g["atac_p_val_adj"] / g["atac_p_val"]).replace([np.inf, -np.inf], np.nan)
        r = r[g["atac_p_val_adj"] < 0.9].dropna()
        uni.append({"tissue": tis, "contrast": con, "celltype": ct,
                    "n_peaks_da": len(g),
                    "n_tests_est": float(np.median(r)) if len(r) else np.nan})
    universe = pd.DataFrame(uni)
    say("\nTested peak universe per tissue (Bonferroni back-calculation):")
    for tis, g in universe.groupby("tissue"):
        v = g["n_tests_est"].dropna()
        if len(v):
            say("  %-6s median %.0f, range %.0f to %.0f over %d cell types"
                % (tis, np.median(v), v.min(), v.max(), len(v)))
            if v.max() - v.min() > 10:
                warn("%s: peak universe is not constant across cell types; "
                     "per-cell-type denominators are needed for proportions" % tis)
else:
    universe = pd.DataFrame()


# %%
# -----------------------------------------------------------------------------
# Cell 5: ALLGENES peak -> gene map
#
# Ragged whitespace: first field is the peak, remaining fields are gene symbols.
# One peak can map to several genes, so any count of DA peaks must deduplicate
# on peak_id rather than counting rows here.
# -----------------------------------------------------------------------------

def read_ragged(path):
    """Yield (key, [partners]) for each non-empty line."""
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                f = line.split()
                if f:
                    yield f[0], f[1:]
    except Exception as exc:
        warn("could not read %s (%s)" % (os.path.relpath(path, ROOT), exc))


say("\n" + "-" * 78)
say("Loading %d ALLGENES peak->gene maps" % len(src_asg_peak))
say("-" * 78)

rows = []
for (tis, con, ct), path in sorted(src_asg_peak.items()):
    if not keep(tis, con):
        continue
    for peak_tok, genes in read_ragged(path):
        p = parse_peak(peak_tok)
        if p is None:
            continue
        if not genes:
            rows.append({"tissue": tis, "contrast": con, "celltype": ct,
                         "peak_id": p[0], "gene": None, "n_genes_for_peak": 0})
            continue
        for g in genes:
            rows.append({"tissue": tis, "contrast": con, "celltype": ct,
                         "peak_id": p[0], "gene": g, "n_genes_for_peak": len(genes)})

peak_gene = pd.DataFrame(rows)
say("\npeak_gene_map: %d peak-gene pairs, %d unique peaks, %d unique genes"
    % (len(peak_gene),
       peak_gene["peak_id"].nunique() if len(peak_gene) else 0,
       peak_gene["gene"].nunique() if len(peak_gene) else 0))

if len(peak_gene):
    multi = peak_gene[peak_gene["n_genes_for_peak"] > 1]
    say("peaks assigned to more than one gene: %d of %d (%.1f%%)"
        % (multi["peak_id"].nunique(), peak_gene["peak_id"].nunique(),
           100 * multi["peak_id"].nunique() / max(peak_gene["peak_id"].nunique(), 1)))


# %%
# -----------------------------------------------------------------------------
# Cell 6: TSS distance from the wide gene-peak matrix
#
# The wide matrix is the only place the signed distance lives. Keyed on
# (tissue, contrast, celltype, gene, peak_id) so it joins onto peak_gene.
# -----------------------------------------------------------------------------

say("\n" + "-" * 78)
say("Loading wide gene-peak matrices for TSS distance")
say("-" * 78)

rows = []
for (tis, con), paths in sorted(src_gpm_wide.items()):
    if not keep(tis, con):
        continue
    # prefer the ALLGENES copy over the root-level duplicate
    path = sorted(paths, key=lambda p: ("ALLGENES" not in p, len(p)))[0]
    say("  %s / %s: %s" % (tis, con, os.path.relpath(path, ROOT)))
    try:
        wide = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as exc:
        warn("could not read %s (%s)" % (os.path.relpath(path, ROOT), exc))
        continue

    gene_col = wide.columns[0]
    for col in wide.columns[1:]:
        m = DA_RE.match(str(col))
        if not m:
            continue
        ct = norm_celltype(m.group("celltype"))
        sub = wide[[gene_col, col]].dropna()
        sub = sub[sub[col] != "NP"]
        for gene, cell in sub.itertuples(index=False):
            for tok in str(cell).split(";"):
                p = parse_peak(tok)
                if p is None:
                    continue
                rows.append({
                    "tissue": tis, "contrast": con, "celltype": ct,
                    "gene": gene, "peak_id": p[0],
                    "distance_to_tss": p[4],
                })

dist = pd.DataFrame(rows)
if len(dist):
    dist = dist.drop_duplicates(
        subset=["tissue", "contrast", "celltype", "gene", "peak_id"])
    say("\ndistance records: %d" % len(dist))
    d = dist["distance_to_tss"].dropna()
    if len(d):
        say("observed TSS distance window: %d to %d  (vendor convention: "
            "negative is upstream)" % (d.min(), d.max()))
        if d.min() > OBSERVED_WINDOW[0] * 1.5 or d.max() < OBSERVED_WINDOW[1] * 0.5:
            warn("distance window differs from the previously observed "
                 "-20000 to +2000; re-check the annotation settings")
        warn("the annotation window is roughly 20 kb upstream to 2 kb downstream, "
             "so enhancers outside that range are absent from these files by "
             "construction, not by absence of signal")
    dist["region"] = dist["distance_to_tss"].map(
        lambda x: classify_region(x, PROMOTER_BP, PROXIMAL_BP))
    say("\nregion breakdown:")
    for k, v in dist["region"].value_counts(dropna=False).items():
        say("  %-20s %7d" % (k, v))
else:
    warn("no TSS distances recovered; downstream region classification will be empty")

peak_gene = (peak_gene.merge(
    dist, on=["tissue", "contrast", "celltype", "gene", "peak_id"], how="left")
    if len(peak_gene) and len(dist) else peak_gene)

if len(peak_gene) and "distance_to_tss" in peak_gene.columns:
    n_nodist = peak_gene["gene"].notna().sum() - peak_gene["distance_to_tss"].notna().sum()
    if n_nodist > 0:
        warn("%d peak-gene pairs from ASSIGNED have no distance in the wide matrix; "
             "the two were assembled by different steps so the sets are not identical"
             % n_nodist)


# %%
# -----------------------------------------------------------------------------
# Cell 7: per-cell-type RNA statistics
#
# logFC here is EXPRESSION. Column names are prefixed rna_ so nothing
# downstream can confuse it with atac_log2FC.
# -----------------------------------------------------------------------------

say("\n" + "-" * 78)
say("Loading %d per-celltype RNA tables" % len(src_gpm_long))
say("-" * 78)

RNA_RENAME = {
    "logFC": "rna_logFC", "logCPM": "rna_logCPM", "F": "rna_F",
    "p_val": "rna_p_val", "p_adj.loc": "rna_p_adj_loc",
    "p_adj.glb": "rna_p_adj_glb",
}

rows = []
for (tis, con, ct), path in sorted(src_gpm_long.items()):
    if not keep(tis, con):
        continue
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as exc:
        warn("could not read %s (%s)" % (os.path.relpath(path, ROOT), exc))
        continue

    df = df.rename(columns=RNA_RENAME)
    gene_col = df.columns[0]
    peak_col = next((c for c in df.columns if str(c).startswith("DA_")), None)

    out = pd.DataFrame({
        "tissue": tis, "contrast": con, "celltype": ct,
        "gene": df[gene_col].astype(str),
    })
    for c in RNA_RENAME.values():
        out[c] = df[c] if c in df.columns else np.nan
    out["rna_tested"] = True
    out["peak_string"] = df[peak_col] if peak_col else None
    out["has_linked_peak"] = (out["peak_string"].notna()
                              & (out["peak_string"] != "NP")) if peak_col else False
    if "contrast" in df.columns:
        out["vendor_contrast_label"] = df["contrast"]
    rows.append(out)

genes_rna = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

if len(genes_rna):
    say("\ngenes_rna: %d rows, %d unique genes, %d groups"
        % (len(genes_rna), genes_rna["gene"].nunique(),
           genes_rna.groupby(["tissue", "contrast", "celltype"]).ngroups))

    # the tested universe should be contrast-independent within a cell type
    sizes = (genes_rna.groupby(["tissue", "celltype", "contrast"])
             .size().unstack(fill_value=0))
    inconsistent = sizes[sizes.replace(0, np.nan).nunique(axis=1) > 1]
    if len(inconsistent):
        warn("tested gene universe differs across contrasts for %d cell types; "
             "it should be contrast-independent, so check these: %s"
             % (len(inconsistent), list(inconsistent.index)[:6]))
    else:
        say("tested gene universe is consistent across contrasts within each "
            "cell type, as expected for a post-filterByExpr set")

    # which cell types have ATAC DA but no RNA table at all
    have_rna = set(map(tuple, genes_rna[["tissue", "contrast", "celltype"]]
                       .drop_duplicates().values))
    have_atac = set(k for k in src_da if keep(k[0], k[1]))
    for k in sorted(have_atac - have_rna):
        warn("%s / %s / %s has ATAC DA results but no RNA table, so it can "
             "appear in peak-only and motif figures but not in the DEG "
             "integration" % k)


# %%
# -----------------------------------------------------------------------------
# Cell 8: integrate
#
# Outer joins throughout. Three kinds of row survive and each is labeled:
#   both      a tested gene with a linked DA peak
#   rna_only  a tested gene with no linked peak
#   atac_only a DA peak whose gene is absent from the tested RNA universe
# -----------------------------------------------------------------------------

say("\n" + "-" * 78)
say("Integrating")
say("-" * 78)

KEY = ["tissue", "contrast", "celltype"]

if len(peak_gene) and len(peaks_da):
    integrated = peak_gene.merge(peaks_da, on=KEY + ["peak_id"], how="outer")
else:
    integrated = peak_gene if len(peak_gene) else peaks_da

if len(integrated) and len(genes_rna):
    integrated = integrated.merge(
        genes_rna.drop(columns=["peak_string"], errors="ignore"),
        on=KEY + ["gene"], how="outer")

if len(integrated):
    integrated["rna_tested"] = integrated.get("rna_tested", pd.Series(dtype=object)).fillna(False)
    has_atac = integrated["atac_log2FC"].notna() if "atac_log2FC" in integrated else False
    integrated["evidence"] = np.where(
        integrated["rna_tested"].astype(bool) & has_atac, "both",
        np.where(integrated["rna_tested"].astype(bool), "rna_only", "atac_only"))

    # concordance, deliberately left unoriented until ident.1 is confirmed
    if "rna_logFC" in integrated and "atac_log2FC" in integrated:
        both = integrated["evidence"] == "both"
        integrated["same_sign_raw"] = np.where(
            both & integrated["rna_logFC"].notna() & integrated["atac_log2FC"].notna(),
            np.sign(integrated["rna_logFC"]) == np.sign(integrated["atac_log2FC"]),
            None)

    integrated["atac_ident1_group"] = ATAC_IDENT1_GROUP
    if ATAC_IDENT1_GROUP is None:
        warn("ATAC_IDENT1_GROUP is unset, so the sign of atac_log2FC cannot be "
             "read as opened or closed. same_sign_raw is a raw sign agreement "
             "only and must not be labeled concordant or discordant in a figure "
             "until Lakmini confirms which group was ident.1")

    say("\nevidence breakdown:")
    for k, v in integrated["evidence"].value_counts().items():
        say("  %-10s %8d" % (k, v))

    say("\nper-group summary:")
    summary = (integrated.groupby(KEY)
               .agg(n_rows=("evidence", "size"),
                    n_peaks=("peak_id", "nunique"),
                    n_genes=("gene", "nunique"),
                    n_both=("evidence", lambda s: (s == "both").sum()),
                    n_rna_only=("evidence", lambda s: (s == "rna_only").sum()),
                    n_atac_only=("evidence", lambda s: (s == "atac_only").sum()))
               .reset_index())
    say(summary.to_string(index=False))
else:
    summary = pd.DataFrame()
    warn("nothing to integrate")


# %%
# -----------------------------------------------------------------------------
# Cell 9: validation
# -----------------------------------------------------------------------------

say("\n" + "-" * 78)
say("Validation")
say("-" * 78)

validation = []

# (a) DA peaks that never appear in the ASSIGNED map
if len(peaks_da) and len(peak_gene):
    for k, g in peaks_da.groupby(KEY):
        mapped = set(peak_gene[
            (peak_gene["tissue"] == k[0]) & (peak_gene["contrast"] == k[1])
            & (peak_gene["celltype"] == k[2])]["peak_id"])
        da_set = set(g["peak_id"])
        validation.append({
            "tissue": k[0], "contrast": k[1], "celltype": k[2],
            "n_da_peaks": len(da_set),
            "n_assigned_peaks": len(mapped),
            "n_da_not_assigned": len(da_set - mapped),
            "n_assigned_not_da": len(mapped - da_set),
        })
    val = pd.DataFrame(validation)
    say("\nDA peaks versus ASSIGNED map coverage:")
    say(val.to_string(index=False))
    tot_unassigned = val["n_da_not_assigned"].sum()
    if tot_unassigned:
        say("\n%d DA peaks across all groups have no gene within the annotation "
            "window. That is expected for intergenic peaks and they are kept as "
            "atac_only rows." % tot_unassigned)
    extra = val[val["n_assigned_not_da"] > 0]
    if len(extra):
        warn("%d groups have peaks in the ASSIGNED map that are not in the DA "
             "table; the two were generated at different times so confirm which "
             "peak set is authoritative" % len(extra))
else:
    val = pd.DataFrame()

# (b) reciprocal agreement between the two ASSIGNED orientations
recip = pd.DataFrame()
if VALIDATE_RECIPROCAL and src_asg_gene:
    say("\nChecking all-peak against all-gene reciprocity ...")
    rows = []
    for key, gpath in sorted(src_asg_gene.items()):
        if not keep(key[0], key[1]) or key not in src_asg_peak:
            continue
        from_gene = set()
        for gene, peaks in read_ragged(gpath):
            for tok in peaks:
                p = parse_peak(tok)
                if p:
                    from_gene.add((p[0], gene))
        from_peak = set()
        for peak_tok, genes in read_ragged(src_asg_peak[key]):
            p = parse_peak(peak_tok)
            if p:
                for g in genes:
                    from_peak.add((p[0], g))
        rows.append({
            "tissue": key[0], "contrast": key[1], "celltype": key[2],
            "n_from_gene": len(from_gene), "n_from_peak": len(from_peak),
            "n_shared": len(from_gene & from_peak),
            "n_gene_only": len(from_gene - from_peak),
            "n_peak_only": len(from_peak - from_gene),
        })
    recip = pd.DataFrame(rows)
    if len(recip):
        say(recip.to_string(index=False))
        bad = recip[(recip["n_gene_only"] > 0) | (recip["n_peak_only"] > 0)]
        if len(bad):
            warn("%d of %d groups show asymmetry between the two ASSIGNED "
                 "orientations; all-peak keeps only the CLOSEST gene per the "
                 "vendor note, so all-gene is the more permissive map"
                 % (len(bad), len(recip)))
        else:
            say("\nthe two orientations agree exactly in every group")

# (c) floor rows per group
floor_tab = pd.DataFrame()
if len(peaks_da):
    floor_tab = (peaks_da.groupby(KEY)
                 .agg(n_peaks=("peak_id", "size"),
                      n_floor=("floor_flag", "sum"),
                      frac_floor=("floor_flag", "mean"),
                      lfc_min_real=("atac_log2FC",
                                    lambda s: s[~peaks_da.loc[s.index, "floor_flag"]].min()),
                      lfc_max_real=("atac_log2FC",
                                    lambda s: s[~peaks_da.loc[s.index, "floor_flag"]].max()))
                 .reset_index())
    say("\nPseudocount-floor rows per group (MIN_PCT_FLOOR=%.3f, "
        "MAX_PLAUSIBLE_ABS_LOG2FC=%.1f):" % (MIN_PCT_FLOOR, MAX_PLAUSIBLE_ABS_LOG2FC))
    say(floor_tab.to_string(index=False))

    heavy = floor_tab[floor_tab["frac_floor"] > 0.2]
    for r in heavy.itertuples():
        warn("%s / %s / %s is %.0f%% floor rows; a volcano or ranked bar plot "
             "here will be dominated by pseudocount artifacts"
             % (r.tissue, r.contrast, r.celltype, 100 * r.frac_floor))

    # single-direction groups, which cannot support an up-versus-down comparison
    for k, g in peaks_da[~peaks_da["floor_flag"]].groupby(KEY):
        s = g["atac_log2FC"].dropna()
        if len(s) and (s > 0).all():
            warn("%s / %s / %s: every non-floor peak moves one way (all positive, "
                 "n=%d); no up-versus-down split or paired motif comparison is "
                 "possible here" % (k + (len(s),)))
        elif len(s) and (s < 0).all():
            warn("%s / %s / %s: every non-floor peak moves one way (all negative, "
                 "n=%d); no up-versus-down split or paired motif comparison is "
                 "possible here" % (k + (len(s),)))


# %%
# -----------------------------------------------------------------------------
# Cell 10: write
# -----------------------------------------------------------------------------

def dump(df, name):
    if df is None or not len(df):
        say("  skipped %s (empty)" % name)
        return
    try:
        df.to_csv(os.path.join(OUT_DIR, name), index=False, lineterminator="\n")
        say("  wrote %-34s %8d rows" % (name, len(df)))
    except Exception as exc:
        warn("could not write %s (%s)" % (name, exc))


say("\n" + "=" * 78)
say("Writing to %s" % OUT_DIR)
say("=" * 78)

dump(peaks_da, "peaks_da.csv")
dump(peak_gene, "peak_gene_map.csv")
dump(genes_rna, "genes_rna.csv")
dump(integrated, "integrated.csv")
dump(summary, "group_summary.csv")
dump(universe, "peak_universe_estimate.csv")
dump(val, "validation_da_vs_assigned.csv")
dump(recip, "validation_reciprocal.csv")
dump(floor_tab, "floor_summary.csv")
dump(pd.DataFrame({"warning": WARNINGS}), "warnings.csv")

meta = {
    "generated": pd.Timestamp.now().isoformat(),
    "root": ROOT,
    "tissues_filter": TISSUES or "all",
    "contrasts_filter": CONTRASTS or "all",
    "min_pct_floor": MIN_PCT_FLOOR,
    "max_plausible_abs_log2fc": MAX_PLAUSIBLE_ABS_LOG2FC,
    "promoter_bp": PROMOTER_BP,
    "proximal_bp": PROXIMAL_BP,
    "atac_ident1_group": ATAC_IDENT1_GROUP,
    "atac_direction_resolved": ATAC_IDENT1_GROUP is not None,
    "n_warnings": len(WARNINGS),
    "note": ("atac_log2FC is accessibility from Signac FindMarkers. rna_logFC is "
             "expression from muscat/edgeR. They are different measurements and "
             "must never share an axis label."),
}
try:
    with open(os.path.join(OUT_DIR, "run_metadata.json"), "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    say("  wrote run_metadata.json")
except Exception as exc:
    warn("could not write run_metadata.json (%s)" % exc)

say("\n%d warnings. Read warnings.csv before using any of this." % len(WARNINGS))
say("Done.")
