# =============================================================================
# multiome_common.py   (revision 2)
#
# Shared helpers for the THC/CBD rhesus multiome collaboration. Imported by
# MULTIOME_00, MULTIOME_01, ARC_00 and ARC_01 so a fix to a normalizer fixes
# every script at once.
#
# -----------------------------------------------------------------------------
# WHAT CHANGED IN REVISION 2
# -----------------------------------------------------------------------------
#
# 1. THE GABAergic QUESTION IS RESOLVED. Two independent sources now agree. The
#    loom Clusters attribute carries 18 brain labels including "GABAergic
#    ADARB2+ Neurons", "GABAergic CCK+ Neurons" and "GABAergic SST+/NPY+
#    Neurons", and NO glutamatergic cluster at all; ATAC_all_clusters.xlsx maps
#    its "GABAergic CCK+ Neurons" sheet to DA_THCCBDvsVEH_GlutamtergicNeurons.txt.
#    So the cluster every DA filename calls Glutamtergic is CCK+ GABAergic.
#
#    The TOKEN stays `Glutamatergic_FILENAME` because it is the join key across
#    every tidy table MULTIOME_01 already produced. The corrected biology lives
#    in DISPLAY_LABELS, which is what figures use. Renaming the token would
#    break those joins for no benefit.
#
# 2. LOOM LABEL ALIASES. The loom uses the vendor's final publication names
#    while the DA filenames use their original cluster names. Both now
#    normalize to the same token.
#
# 3. CLUSTER LABELS CONTAIN LITERAL NEWLINES. The loom holds "Smooth\nmuscle
#    cells", "Endothelial\nCells-1", "Entero\ncyte-3",
#    "Entero\nchromaffin\ncells" and "Transient-\namplifying cells". Those
#    corrupt CSV output and break joins. clean_label() collapses whitespace for
#    display; norm_celltype() strips all non-alphanumerics, which matters
#    because the newline sits MID-WORD: "Entero\ncyte-3" must become
#    "enterocyte3", not "Entero cyte-3".
#
# 4. CELL IDENTIFIER PARSING. parse_cell_id() handles the real loom format,
#    BARCODE-1_LIBRARYINDEX, which is Seurat merge() appending a list position
#    to a 10x barcode that already carries its own -1.
#
# 5. CLUSTERS WITH NO DA RESULTS are recorded, so a missing cell type is a
#    known fact rather than a surprise mid-figure.
#
# Jake Lehle, 2026-10
# =============================================================================

import os
import re

__all__ = [
    "CONTRAST_ALIASES", "CELLTYPE_ALIASES", "DISPLAY_LABELS",
    "BRAIN_LABEL_WARNINGS", "NO_DA_CLUSTERS", "TREATMENT_ARMS",
    "norm_contrast", "norm_celltype", "clean_label", "display_label",
    "infer_tissue", "parse_cell_id",
    "parse_peak", "classify_region", "peek_gpm_header",
    "FASTQ_RE", "DA_RE", "ASSIGNED_RE", "OBSERVED_WINDOW",
    "VENDOR_TSS_CONVENTION",
]


# -----------------------------------------------------------------------------
# Resolved facts, recorded here so no script has to rediscover them
# -----------------------------------------------------------------------------

# From the loom Stimulation attribute cross-tabulated against Samples, checked
# against the arms in Sample List.xlsx. Arm A is THC/CBD, Arm B is VEH.
TREATMENT_ARMS = {
    "SIV/cART/THCCBD":  {"40872", "40874", "40877", "40878"},   # brain spelling
    "SIV/cART/THC:CBD": {"40872", "40874", "40877", "40878"},   # colon spelling
    "SIV/cART/VEH":     {"40873", "40875", "40876", "40879"},
    "Pre-infection":    {"40874", "40875", "40877", "40878"},   # colon only
}

# SETUP_03_validate_gtf.py: 100.00% of the vendor's 43,128 unique (peak, gene,
# distance) triples reproduce exactly under this convention, against Ensembl
# release 97 Mmul_8.0.1. Measured, not assumed. Use it for every region
# classification so our numbers stay comparable to theirs.
VENDOR_TSS_CONVENTION = {
    "tss_source": "gene_level",
    "peak_anchor": "midpoint",
    "ensembl_release": 97,
    "assembly": "Mmul_8.0.1",
    "sign": "negative is upstream",
}

# Loom clusters with no DA results. Brain: all small enough that a test was not
# worth running. Colon: Transient-amplifying at 489 cells is NOT small, so that
# omission is worth asking about.
NO_DA_CLUSTERS = {
    "Brain": {
        "Microglia_2": 124, "Microglia_3": 70, "Ependymal": 85,
        "Vascular": 72, "D2M": 37, "Unknown_brain_extra": 21,
    },
    "Colon": {
        "Transient_amplifying": 489,
    },
}


# -----------------------------------------------------------------------------
# Token normalizers
# -----------------------------------------------------------------------------

CONTRAST_ALIASES = {
    "thccbdvspre": "THCCBD_vs_Pre",
    "thccbdvspreinf": "THCCBD_vs_Pre",
    "vehvspre": "VEH_vs_Pre",
    "vehvspreinf": "VEH_vs_Pre",
    "vehvsthccbd": "VEH_vs_THCCBD",
    "vehvsthccbb": "VEH_vs_THCCBD",   # CBB typo in ATAC_DifAccePerGroup_Colon
    "thccbdvsveh": "THCCBD_vs_VEH",
}

CELLTYPE_ALIASES = {
    # ---------------- colon, DA filenames ----------------
    "bcells": "B_cells",
    "bcell": "B_cells",
    "endothelialcells1": "Endothelial_1",
    "endothelialcells2": "Endothelial_2",
    "entericgliacells": "Enteric_glia",
    "enterochromaffincells": "Enterochromaffin",
    "enterocyte1": "Enterocyte_1",
    "enterocyte2": "Enterocyte_2",
    "enterocyte3": "Enterocyte_3",
    "entercyte1": "Enterocyte_1",
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
    # colon loom labels not present in the DA set. The squash normalizer
    # already folds the newline-containing spellings onto the entries above,
    # since it drops all non-alphanumerics.
    "transientamplifyingcells": "Transient_amplifying",

    # ---------------- brain, DA filenames ----------------
    "astrocytes": "Astrocytes",
    "d1mm": "D1M_M",
    "d1ms": "D1M_S",
    "d2mm": "D2M_M",
    "d2ms": "D2M_S",
    "gabaergicneurons1": "GABAergic_1",
    "gabaergicneurons2": "GABAergic_2",
    "glutamtergicneurons": "Glutamatergic_FILENAME",
    "glutamatergicneurons": "Glutamatergic_FILENAME",
    "interneuron": "Interneuron",
    "microglia1": "Microglia_1",
    "microglia2": "Microglia_2",
    "microglia3": "Microglia_3",
    "oligodendrocytes": "Oligodendrocytes",
    "opc": "OPC",

    # ---------------- brain, loom labels ----------------
    # The loom carries the vendor's FINAL names; the DA filenames carry their
    # original ones. Both fold onto the same token so a join works either way.
    "gabaergicsstnpyneurons": "GABAergic_1",
    "gabaergicadarb2neurons": "GABAergic_2",
    "gabaergiccckneurons": "Glutamatergic_FILENAME",   # the resolved one
    "tac3interneurons": "Interneuron",
    "opcs": "OPC",
    "ependymalcells": "Ependymal",
    "vascularcells": "Vascular",
    "d2m": "D2M",
}

# What figures should print. Tokens are join keys and stay stable; these are
# the biology.
DISPLAY_LABELS = {
    # brain
    "Astrocytes": "Astrocytes",
    "D1M_M": "D1M-M",
    "D1M_S": "D1M-S",
    "D2M_M": "D2M-M",
    "D2M_S": "D2M-S",
    "D2M": "D2M",
    "GABAergic_1": "GABAergic SST+/NPY+ Neurons",
    "GABAergic_2": "GABAergic ADARB2+ Neurons",
    # Resolved: the DA filenames call this Glutamtergic, but the loom and
    # ATAC_all_clusters.xlsx both call it CCK+ GABAergic, and the loom has no
    # glutamatergic cluster at all. Striatum is overwhelmingly GABAergic.
    "Glutamatergic_FILENAME": "GABAergic CCK+ Neurons",
    "Interneuron": "TAC3+ Interneurons",
    "Microglia_1": "Microglia-1",
    "Microglia_2": "Microglia-2",
    "Microglia_3": "Microglia-3",
    "Oligodendrocytes": "Oligodendrocytes",
    "OPC": "OPC",
    "Ependymal": "Ependymal cells",
    "Vascular": "Vascular cells",
    # colon
    "B_cells": "B cells",
    "T_cells": "T cells",
    "NK_cells": "NK cells",
    "Plasma_cells": "Plasma cells",
    "Macrophage": "Macrophage",
    "Endothelial_1": "Endothelial Cells-1",
    "Endothelial_2": "Endothelial Cells-2",
    "Enteric_glia": "Enteric glia cells",
    "Enterochromaffin": "Enterochromaffin cells",
    "Enterocyte_1": "Enterocyte-1",
    "Enterocyte_2": "Enterocyte-2",
    "Enterocyte_3": "Enterocyte-3",
    "Fibroblast_S1a": "Fibroblast S1a",
    "Fibroblast_S1b": "Fibroblast S1b",
    "Fibroblast_S1c": "Fibroblast S1c",
    "Fibroblast_S2a": "Fibroblast S2a",
    "Fibroblast_S2b": "Fibroblast S2b",
    "Smooth_muscle": "Smooth muscle cells",
    "Transient_amplifying": "Transient-amplifying cells",
    "Unknown": "Unknown",
}

# Kept for the record: the discrepancy between DA filename tokens and the
# vendor's final sheet names, now corroborated by the loom.
BRAIN_LABEL_WARNINGS = {
    "GABAergic_1": "GABAergic SST+ NPY+",
    "GABAergic_2": "GABAergic ADARB2+ neurons",
    "Glutamatergic_FILENAME": "GABAergic CCK+ Neurons",
    "Interneuron": "TAC3+ Interneurons",
    "Microglia_1": "Microglia",
    "OPC": "OPCs",
}


def clean_label(s):
    """
    Collapse the literal newlines and runs of whitespace the loom carries in
    its cluster labels, for display and for safe CSV writing.

    Note the newline sits mid-word, so this gives "Entero cyte-3" rather than
    "Enterocyte-3". Use norm_celltype() for joining, which drops all
    non-alphanumerics and therefore recovers "enterocyte3".
    """
    return re.sub(r"\s+", " ", str(s).replace("\n", " ").replace("\r", " ")).strip()


def _squash(token):
    """Lowercase, then drop everything that is not a letter or digit."""
    return re.sub(r"[^a-z0-9]", "", str(token).lower())


def norm_contrast(token):
    return CONTRAST_ALIASES.get(_squash(token), "UNKNOWN:" + str(token))


def norm_celltype(token):
    """
    Fold a cell type label from any source (DA filename, cluster_id column,
    loom Clusters attribute) onto one join token. Robust to the embedded
    newlines, because squashing drops them along with spaces and punctuation.
    """
    return CELLTYPE_ALIASES.get(_squash(token), "UNKNOWN:" + str(token))


def display_label(token):
    """Publication label for a token, falling back to the token itself."""
    return DISPLAY_LABELS.get(token, token)


# -----------------------------------------------------------------------------
# Cell identifiers
# -----------------------------------------------------------------------------

# The real loom format: AAACAGCCACTATGGC-1_1
#   AAACAGCCACTATGGC   16 bp 10x barcode
#   -1                 the per-library suffix 10x itself writes
#   _1                 the list position Seurat's merge() appends
#
# So the LIBRARY INDEX is the trailing integer after the last underscore, and
# that index is what an aggr CSV ordering has to reproduce.
CELL_ID_RE = re.compile(
    r"^(?P<barcode>[ACGTN]+)(?:-(?P<tenx>\d+))?(?:_(?P<lib>\d+))?$"
)


def parse_cell_id(cell_id):
    """
    'AAACAGCCACTATGGC-1_1' -> ('AAACAGCCACTATGGC', 1, 1)
    'AAACAGCCACTATGGC-1'   -> ('AAACAGCCACTATGGC', 1, None)
    'AAACAGCCACTATGGC'     -> ('AAACAGCCACTATGGC', None, None)

    Returns (barcode, tenx_suffix, library_index) or (None, None, None).
    """
    m = CELL_ID_RE.match(str(cell_id).strip())
    if not m:
        return (None, None, None)
    bc = m.group("barcode")
    tenx = int(m.group("tenx")) if m.group("tenx") else None
    lib = int(m.group("lib")) if m.group("lib") else None
    return (bc, tenx, lib)


# -----------------------------------------------------------------------------
# Tissue inference
#
# The FASTQ tree sits under "Raw Data_Colon_Brain", whose name contains BOTH
# tissue words, so any substring test on the full path is wrong. That bug
# labelled all 12 colon libraries as Brain. Resolution order: the explicit
# "Brain data" / "Colon data" component, then a top-level directory naming
# exactly one tissue, then the sample-name prefix.
# -----------------------------------------------------------------------------

def infer_tissue(parts):
    for p in parts:
        pl = str(p).lower()
        if pl in ("brain data", "colon data"):
            return "Brain" if pl.startswith("brain") else "Colon"

    top = str(parts[0]).lower() if parts else ""
    has_b, has_c = "brain" in top, "colon" in top
    if has_b and not has_c:
        return "Brain"
    if has_c and not has_b:
        return "Colon"

    base = str(parts[-1]) if parts else ""
    if re.match(r"^BG_\d+", base):
        return "Brain"
    if re.match(r"^C\d{5}_[NP][AG]", base):
        return "Colon"
    return None


# -----------------------------------------------------------------------------
# Filename patterns
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

GPM_HEADER_RE = re.compile(
    r"^DA_(?P<contrast>[^_]+)_(?P<celltype>.+)\.txt$", re.IGNORECASE)


def peek_gpm_header(path):
    """
    Read only the header of a gene-peak-matrix-peaks.txt and return
    (shape, contrast, celltypes), where shape is "wide" or "long".

    The file is self-describing: its DA_* column names carry both contrast and
    cell type, which is the only reliable source for brain, where no contrast
    directory exists in the path.
    """
    try:
        with open(path, "r", errors="replace") as fh:
            header = fh.readline().rstrip("\n").rstrip("\r")
    except Exception:
        return (None, None, [])
    if not header:
        return (None, None, [])

    fields = header.split("\t")
    da_fields = [f for f in fields if f.startswith("DA_")]
    if not da_fields:
        return (None, None, [])

    contrasts, celltypes = set(), []
    for f in da_fields:
        m = GPM_HEADER_RE.match(f)
        if not m:
            continue
        contrasts.add(norm_contrast(m.group("contrast")))
        celltypes.append(norm_celltype(m.group("celltype")))

    contrast = contrasts.pop() if len(contrasts) == 1 else None
    shape = "wide" if len(da_fields) >= 2 else "long"
    return (shape, contrast, celltypes)


# -----------------------------------------------------------------------------
# Peak coordinates
#
#     chr1-56999-57879          bare, in DA and ASSIGNED files
#     chr1-56999-57879_-25      with signed TSS distance, in the gene-peak matrix
#
# rsplit throughout so a contig name containing '-' or '_' and a negative
# distance both survive. The delivered data has only chr1-chr20, chrX, chrY,
# but the parser stays general.
# -----------------------------------------------------------------------------

_PEAK_WITH_DIST_RE = re.compile(r"^(?P<peak>.+?)_(?P<dist>-?\d+)$")


def parse_peak(token):
    """
    Return (peak_id, chrom, start, end, distance) or None. peak_id excludes any
    distance suffix so it joins cleanly against the DA and ASSIGNED files.
    """
    token = str(token).strip()
    if not token or token in ("NP", "nan", "NA", "None"):
        return None

    dist = None
    m = _PEAK_WITH_DIST_RE.match(token)
    if m:
        token = m.group("peak")
        dist = int(m.group("dist"))

    bits = token.rsplit("-", 2)
    if len(bits) != 3:
        return None
    chrom, start, end = bits
    try:
        return (token, chrom, int(start), int(end), dist)
    except ValueError:
        return None


# -----------------------------------------------------------------------------
# Region classification
#
# Negative distance is upstream, per the vendor note and confirmed by
# SETUP_03's 100% reproduction. Observed window is about -20,000 to +2,000, so
# anything outside it was never annotated: a null there is not evidence against
# a distal element.
# -----------------------------------------------------------------------------

OBSERVED_WINDOW = (-20000, 2000)


def classify_region(distance, promoter_bp=1000, proximal_bp=10000):
    if distance is None:
        return None
    d = int(distance)
    if abs(d) <= promoter_bp:
        return "promoter"
    if d < 0:
        return "upstream_proximal" if d >= -proximal_bp else "upstream_distal"
    return "downstream"
