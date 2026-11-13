#!/usr/bin/env Rscript
# =============================================================================
# ARC_06_select_celltypes.R  <tissue> [select|full]
#
# Answers "which cell type" for the 22 genes in the collaborators' brain
# multiome report, then hands the surviving gene and cell-type pairs to
# ARC_05 as a panel.
#
#   Rscript ARC_06_select_celltypes.R Brain select   # phases 0-2, seconds, no object
#   Rscript ARC_06_select_celltypes.R Brain full     # all phases, loads the object
#   Rscript ARC_06_select_celltypes.R Colon full
#
# select mode reads only the tidy CSVs and answers the question you actually
# asked: which cell types show these genes moving. It takes seconds and needs
# no memory. Run it first, read deg_matrix_<tissue>.csv, and only then spend
# an hour on full.
#
# full mode adds the peak inventory and the per-animal pseudobulk recompute,
# both of which need seurat_<tissue>_rebuilt.rds.
#
# -----------------------------------------------------------------------------
# THE UNIT OF ANALYSIS IS gene x cell type x contrast
#
# Tissue only selects which object to open. Every test is run inside one cell
# type. Nothing is pooled across cell types, because a gene moving in microglia
# and not in oligodendrocytes is the answer, not noise to average away.
# -----------------------------------------------------------------------------
#
# WHAT WAS VERIFIED AGAINST THE SOURCE SCRIPTS BEFORE THIS WAS WRITTEN
#
# MULTIOME_01_build_tidy.py, OUT_DIR at line 57:
#   /master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS/01_tidy
#
# genes_rna.csv columns (built at MULTIOME_01 cell 7):
#   tissue, contrast, celltype, gene, rna_logFC, rna_logCPM, rna_F,
#   rna_p_val, rna_p_adj_loc, rna_p_adj_glb, rna_tested, peak_string,
#   has_linked_peak, vendor_contrast_label
#
# There is NO column called FDR. rna_p_adj_loc is muscat's p_adj.loc, which is
# already corrected WITHIN each cluster, and is therefore the right one for a
# per-cell-type test. rna_p_adj_glb corrects across all clusters at once and is
# reported alongside but not used to select.
#
# peaks_da.csv columns (MULTIOME_01 cell 4):
#   tissue, contrast, celltype, peak_id, chrom, start, end, width,
#   atac_log2FC, atac_p_val, atac_p_val_adj, pct_1, pct_2, min_pct,
#   floor_flag, atac_direction_raw
#
# Coordinates are already split out, so nothing here regex-parses a peak string.
#
# Contrast tokens come from CONTRAST_ALIASES in multiome_common.py:
#   VEH_vs_THCCBD, THCCBD_vs_VEH, VEH_vs_Pre, THCCBD_vs_Pre
#
# celltype in the tidy tables is a TOKEN (Microglia_1, GABAergic_1). The Seurat
# object carries celltype_display (Microglia-1, GABAergic SST+/NPY+ Neurons).
# DISPLAY_LABELS in multiome_common.py is the only bridge, and several pairs
# cannot be derived by normalizing, so the map is copied here and gated.
#
# Jake Lehle, 2026-11
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Cell 1: configuration
# -----------------------------------------------------------------------------

ANALYSIS_ROOT <- "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
TIDY_DIR      <- file.path(ANALYSIS_ROOT, "01_tidy")
OBJ_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "04_objects")
OUT_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "06_selection")

# THE filter. Deliberately stringent: we want the cell types with a large
# effect in the collaborators' top picks, and an empty result is a question for
# them rather than a reason to loosen. Absolute value, not one-sided, because
# the report's own claim about decreased interferon-stimulated genes would be
# thrown away by a one-sided cut.
LFC_PRIMARY <- 1.0
FDR_PRIMARY <- 0.05

# Reported alongside, never used to select. The sweep moves BOTH axes, because
# the first brain run showed the binding constraint is FDR and not logFC:
# VSIG4 reached logFC 3.50 in Microglia-1 at FDR 0.913, and the best FDR
# anywhere in the panel was 0.472. A sweep that holds FDR at 0.05 and varies
# logFC returns zeros and diagnoses nothing.
LFC_SWEEP <- c(0.25, 0.50, 1.00)
FDR_SWEEP <- c(0.05, 0.10, 0.25)

# Which multiple-testing correction selection runs on.
#   "loc"   rna_p_adj_loc, BH within cluster across EVERY gene tested there
#           (about 7,000 in brain). The vendor's own column.
#   "glb"   rna_p_adj_glb, BH across all clusters at once. Stricter still.
#   "panel" BH computed here, within cluster, across ONLY the panel genes
#           tested in that cluster. Defensible because the 23 genes are a
#           targeted list handed to us by the collaborators rather than a
#           screen, so correcting across 7,000 genes answers a question nobody
#           asked. It is a weaker claim than "loc" and must be labelled as a
#           targeted follow-up wherever it is reported.
# All three are always computed and compared. This only sets which one gates.
# Changed from "loc" to "panel" on 2026-11-05 after the first brain full run.
# Not a loosened threshold. At FDR_BASIS "loc" the best value in the whole
# panel was 0.472, and capturing NRP1, VSIG4 and FOLR2 would have required a
# threshold of 0.92, which is not error control. Their raw p-values are
# 2.8e-4, 0.011 and 0.018; they only looked bad because BH was correcting them
# against roughly 7,000 genes per cluster. "panel" changes the denominator to
# the ~20 genes the collaborators actually asked about, which is the question
# being posed. Read the circularity warning this triggers before reporting it.
FDR_BASIS <- "panel"

# Three names, because one name cannot do three jobs and trying cost a run.
# The source table genes_rna.csv calls these columns rna_p_adj_loc and
# rna_p_adj_glb. This script renames them to fdr_loc and fdr_glb when it builds
# its own tables, and fdr_panel is computed here and exists nowhere upstream.
# So: SRC for reading and validating the source table, GATE for anything
# operating on prov/hits, LABEL for anything a human reads.
FDR_GATE_COL <- switch(FDR_BASIS, loc = "fdr_loc", glb = "fdr_glb",
                       panel = "fdr_panel",
                       stop("FDR_BASIS must be loc, glb or panel"))
FDR_SRC_COL  <- switch(FDR_BASIS, loc = "rna_p_adj_loc",
                       glb = "rna_p_adj_glb", panel = "rna_p_val")
FDR_LABEL    <- switch(FDR_BASIS,
                       loc   = "rna_p_adj_loc (within cluster, all genes)",
                       glb   = "rna_p_adj_glb (across all clusters)",
                       panel = "fdr_panel (within cluster, panel genes only)")

# Contrast to select on, first match wins. Direction is calibrated below rather
# than assumed from the name.
CONTRAST_PRIMARY <- c("VEH_vs_THCCBD", "THCCBD_vs_VEH")

# Gate C. A pseudobulk contrast needs enough cells per animal within the cell
# type, and enough animals per arm clearing that floor. Not hypothetical: in
# colon macrophages at necropsy the VEH animals hold 6, 2, 16 and 64 cells, so
# at a 10-cell floor only 2 of 4 VEH animals qualify.
MIN_CELLS_PER_ANIMAL <- 10
MIN_ANIMALS_PER_ARM  <- 3

# Only for inspecting a known-broken table. Any run with this TRUE is for
# diagnosis, never for reporting.
ALLOW_DUPLICATE_CELLTYPES <- FALSE

# -----------------------------------------------------------------------------
# CORRUPTED CELL TYPE LABELS
#
# MULTIOME_02_duplicate_audit.py proved that two brain directories carry another
# cell type's numbers under their own correct header. MULTIOME_03 plus Jake's
# earlier independent analysis identified which label is wrong:
#
#   .../Oligodendrocytes/gene-peak-matrix-peaks.txt  holds Microglia-1's data
#   .../Interneuron/gene-peak-matrix-peaks.txt       holds Glutamatergic's data
#
# The marker evidence for the microglia call is detection count, not mean
# logCPM: 8 of 10 microglial markers survive filterByExpr in that file against
# 3 of 10 oligodendrocyte markers, and the three present (MBP, PLP1, MAG) are
# the myelin transcripts that appear as ambient contamination in any microglial
# prep from white matter. MOG, MOBP, CNP, CLDN11, ST18, MYRF and ERMN are all
# absent, which cannot happen in real oligodendrocytes.
#
# These are dropped BEFORE the duplicate check, so the gate still protects
# against a NEW duplication appearing in a future delivery.
# -----------------------------------------------------------------------------
CORRUPTED_CELLTYPES <- list(
  Brain = c("Oligodendrocytes", "Interneuron"),
  Colon = character(0)
)

PROMOTER_BP <- 2000
WINDOW_UP   <- 20000
WINDOW_DOWN <- 5000

FONT_SIZE <- 28
FIG_DPI   <- 300

ARM_COLORS <- c("VEH" = "#0072B2", "THC/CBD" = "#D55E00")
HEAT_LOW   <- "#2166AC"
HEAT_MID   <- "#F7F7F7"
HEAT_HIGH  <- "#B2182B"

# -----------------------------------------------------------------------------
# Cell 2: the gene panel
#
# Transcribed from the table in the collaborators' brain multiome report.
# reported_rna is their claim, used ONLY to calibrate the sign convention of
# rna_logFC, never to filter. reported_atac is their accessibility column.
# -----------------------------------------------------------------------------

GENE_PANEL <- data.table(
  gene = c("NRP1", "VSIG4", "FOLR2", "SIGLEC9", "SIGLEC5",
           "MT-ND1", "MT-ND5", "COX4I1", "NDUFB4", "ATP5MC3",
           "ATP6V0A1", "LAMTOR2", "RNASET2",
           "RPL18", "RPL29", "RPS6", "RPS9", "EEF1A1",
           "UBA52", "NEDD8", "PSMB6", "HSP90AB1",
           "ISG15"),
  pathway = c(rep("Immune regulation", 5),
              rep("Mitochondrial OXPHOS", 5),
              rep("Lysosomal function", 3),
              rep("Ribosomal/proteostasis", 9),
              "Interferon response"),
  reported_rna = c("Up", "Up", "Up", "Up", "Up",
                   "Up", "Up", "Up", "Up", "Up",
                   "Up", "Slightly up", "Slightly up",
                   "Up", "Up", "Up", "Up", "Up", "Up", "Up", "Up", "Up",
                   "Down"),
  reported_atac = c("Yes", "Not reported", "Not reported", "Yes", "Yes",
                    rep("Not reported", 5),
                    "Yes", "Yes", "Yes",
                    rep("Not reported", 9),
                    "Not reported"),
  # Mitochondrial genes can never have peaks here: the reference set
  # non_nuclear_contigs = ["chrM"], cellranger-arc excludes chrM from peak
  # calling, and the consensus set came back on 22 contigs (20 autosomes plus
  # X and Y). Also correct biologically, since mtDNA is not nucleosomal.
  atac_possible = c(rep(TRUE, 5), FALSE, FALSE, rep(TRUE, 3),
                    rep(TRUE, 3), rep(TRUE, 9), TRUE),
  atac_impossible_reason = c(rep(NA_character_, 5),
                             "chrM excluded from peak calling; mtDNA not nucleosomal",
                             "chrM excluded from peak calling; mtDNA not nucleosomal",
                             rep(NA_character_, 3),
                             rep(NA_character_, 3), rep(NA_character_, 9), NA_character_),
  # Housekeeping loci whose promoters are open in nearly every cell type. An
  # open promoter there is close to guaranteed and carries little information.
  weak_target = c(rep(FALSE, 5), rep(FALSE, 5), rep(FALSE, 3),
                  TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE,
                  FALSE)
)

# Aliases tried in order against the tidy tables and the annotation.
GENE_ALIASES <- list(
  SIGLEC9  = c("SIGLEC9", "SIGLEG9", "SIGLEC-9"),
  SIGLEC5  = c("SIGLEC5", "SIGLEG5", "SIGLEC-5"),
  NRP1     = c("NRP1", "NRP-1"),
  ISG15    = c("ISG15", "G1P2", "UCRP", "IFI15", "IGI15"),
  `MT-ND1` = c("MT-ND1", "MTND1", "ND1"),
  `MT-ND5` = c("MT-ND5", "MTND5", "ND5"),
  ATP5MC3  = c("ATP5MC3", "ATP5G3"),
  UBA52    = c("UBA52"),
  HSP90AB1 = c("HSP90AB1", "HSP90B")
)

# -----------------------------------------------------------------------------
# Cell 3: celltype token to display label
#
# Copied verbatim from DISPLAY_LABELS in multiome_common.py. It is duplicated
# rather than imported because this is R and that is Python, and duplication is
# acceptable only because Phase 0 below refuses to run when the two disagree
# with what the object actually contains. Same pattern as LIB_ORDER in ARC_04.
# -----------------------------------------------------------------------------

DISPLAY_LABELS <- c(
  # brain
  "Astrocytes"             = "Astrocytes",
  "D1M_M"                  = "D1M-M",
  "D1M_S"                  = "D1M-S",
  "D2M_M"                  = "D2M-M",
  "D2M_S"                  = "D2M-S",
  "D2M"                    = "D2M",
  "GABAergic_1"            = "GABAergic SST+/NPY+ Neurons",
  "GABAergic_2"            = "GABAergic ADARB2+ Neurons",
  "Glutamatergic_FILENAME" = "GABAergic CCK+ Neurons",
  "Interneuron"            = "TAC3+ Interneurons",
  "Microglia_1"            = "Microglia-1",
  "Microglia_2"            = "Microglia-2",
  "Microglia_3"            = "Microglia-3",
  "Oligodendrocytes"       = "Oligodendrocytes",
  "OPC"                    = "OPC",
  "Ependymal"              = "Ependymal cells",
  "Vascular"               = "Vascular cells",
  # colon
  "B_cells"                = "B cells",
  "T_cells"                = "T cells",
  "NK_cells"               = "NK cells",
  "Plasma_cells"           = "Plasma cells",
  "Macrophage"             = "Macrophage",
  "Endothelial_1"          = "Endothelial Cells-1",
  "Endothelial_2"          = "Endothelial Cells-2",
  "Enteric_glia"           = "Enteric glia cells",
  "Enterochromaffin"       = "Enterochromaffin cells",
  "Enterocyte_1"           = "Enterocyte-1",
  "Enterocyte_2"           = "Enterocyte-2",
  "Enterocyte_3"           = "Enterocyte-3",
  "Fibroblast_S1a"         = "Fibroblast S1a",
  "Fibroblast_S1b"         = "Fibroblast S1b",
  "Fibroblast_S1c"         = "Fibroblast S1c",
  "Fibroblast_S2a"         = "Fibroblast S2a",
  "Fibroblast_S2b"         = "Fibroblast S2b",
  "Smooth_muscle"          = "Smooth muscle cells",
  "Transient_amplifying"   = "Transient-amplifying cells",
  "Unknown"                = "Unknown"
)

# Cell types present in the rebuilt Seurat objects, from the ARC_05 logs.
# Compared against the RNA table so the gap is visible: brain RNA covers only
# 12 of the 18 brain clusters, which matches the report's own statement that
# "differential gene expression data [is] available for 12 clusters". A cluster
# with no RNA statistics can never be selected here, and that includes
# Microglia-2 and Microglia-3.
OBJECT_CELLTYPES <- list(
  Brain = c("Oligodendrocytes", "D2M-M", "Microglia-1", "D1M-M", "Astrocytes",
            "OPC", "D1M-S", "GABAergic ADARB2+ Neurons", "GABAergic CCK+ Neurons",
            "D2M-S", "GABAergic SST+/NPY+ Neurons", "TAC3+ Interneurons",
            "Microglia-2", "Ependymal cells", "Vascular cells", "Microglia-3",
            "D2M", "Unknown"),
  Colon = c("Plasma cells", "T cells", "NK cells", "Fibroblast S1a",
            "Enterocyte-1", "B cells", "Fibroblast S1c", "Enterocyte-2",
            "Smooth muscle cells", "Endothelial Cells-1", "Macrophage",
            "Fibroblast S1b", "Transient-amplifying cells", "Fibroblast S2a",
            "Enterocyte-3", "Endothelial Cells-2", "Unknown",
            "Enteric glia cells", "Fibroblast S2b", "Enterochromaffin cells")
)

# Animal to arm, same table ARC_05 uses, cross-checked there against the loom.
ANIMAL_ARM <- c("40872" = "THC/CBD", "40874" = "THC/CBD",
                "40877" = "THC/CBD", "40878" = "THC/CBD",
                "40873" = "VEH",     "40875" = "VEH",
                "40876" = "VEH",     "40879" = "VEH")

# -----------------------------------------------------------------------------
# Cell 4: plumbing
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: ARC_06_select_celltypes.R <Brain|Colon> [select|full]")
TISSUE <- args[1]
if (!TISSUE %in% c("Brain", "Colon")) stop("tissue must be Brain or Colon")
MODE <- if (length(args) > 1) args[2] else "select"
if (!MODE %in% c("select", "full")) stop("mode must be select or full")

say  <- function(...) cat(paste0(..., "\n"))
hdr  <- function(x) { say(""); say(strrep("=", 78)); say(x); say(strrep("=", 78)) }
subh <- function(x) { say(""); say(paste0("--- ", x, " ---")) }
ok   <- function(x) say(paste0("  [ok]    ", x))
wrn  <- function(x) say(paste0("  [warn]  ", x))

# Same guard as ARC_04 and ARC_05: a helper named sub() shadows base::sub and
# the failure lands deep in a later phase.
for (fn in c("say", "hdr", "subh", "ok", "wrn")) {
  if (exists(fn, envir = baseenv(), inherits = FALSE)) {
    stop("helper '", fn, "' collides with base::", fn, ". Rename it.")
  }
}

OUT_DIR <- file.path(OUT_ROOT, TISSUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

hdr(paste0("ARC_06_select_celltypes : ", TISSUE, "  [", MODE, "]"))
say(paste0("host:     ", Sys.info()[["nodename"]]))
say(paste0("started:  ", format(t_start, "%Y-%m-%d %H:%M:%S")))
say(paste0("tidy:     ", TIDY_DIR))
say(paste0("out:      ", OUT_DIR))
say(paste0("filter:   |rna_logFC| > ", LFC_PRIMARY, " and ", FDR_LABEL,
           " < ", FDR_PRIMARY))
say(paste0("panel:    ", nrow(GENE_PANEL), " genes, ",
           sum(!GENE_PANEL$atac_possible), " of which can never have ATAC peaks"))

dump_csv <- function(dt, name) {
  p <- file.path(OUT_DIR, name)
  fwrite(dt, p)
  say(sprintf("    wrote %s  (%d rows)", name, nrow(dt)))
  invisible(p)
}

# -----------------------------------------------------------------------------
# Cell 5: Phase 0, load and validate
# -----------------------------------------------------------------------------

hdr("Phase 0: load the tidy tables and validate the joins")

rna_path <- file.path(TIDY_DIR, "genes_rna.csv")
da_path  <- file.path(TIDY_DIR, "peaks_da.csv")
if (!file.exists(rna_path)) {
  stop("not found: ", rna_path,
       "\nRun MULTIOME_01_build_tidy.py first. Its OUT_DIR is ", TIDY_DIR)
}
rna <- fread(rna_path)
say(sprintf("  genes_rna.csv: %s rows", format(nrow(rna), big.mark = ",")))

req_rna <- unique(c("tissue", "contrast", "celltype", "gene", "rna_logFC",
                    "rna_p_adj_loc", "rna_p_val", FDR_SRC_COL))
miss <- setdiff(req_rna, names(rna))
if (length(miss)) {
  say("")
  say("  columns present:")
  say(paste0("    ", paste(names(rna), collapse = ", ")))
  stop("genes_rna.csv is missing: ", paste(miss, collapse = ", "),
       "\nThe column layout changed since MULTIOME_01 was written.")
}
ok(paste0("all required columns present, gating on ", FDR_LABEL))

da <- NULL
if (file.exists(da_path)) {
  da <- fread(da_path)
  req_da <- c("tissue", "contrast", "celltype", "peak_id", "chrom", "start", "end")
  if (length(setdiff(req_da, names(da)))) {
    wrn("peaks_da.csv is missing expected columns; vendor DA annotation disabled")
    da <- NULL
  } else {
    ok(sprintf("peaks_da.csv: %s rows, coordinates already split out",
               format(nrow(da), big.mark = ",")))
  }
} else {
  wrn(paste0("no peaks_da.csv at ", da_path, "; vendor DA annotation disabled"))
}

rna_t <- rna[tissue == TISSUE]
if (!nrow(rna_t)) {
  say("")
  say("  tissues present: ", paste(sort(unique(rna$tissue)), collapse = ", "))
  stop("no RNA rows for tissue '", TISSUE, "'")
}

subh("contrasts present for this tissue")
ct_tab <- rna_t[, .(rows = .N, celltypes = uniqueN(celltype),
                    genes = uniqueN(gene)), by = contrast][order(-rows)]
print(ct_tab)

CONTRAST <- intersect(CONTRAST_PRIMARY, unique(rna_t$contrast))
if (!length(CONTRAST)) {
  stop("none of ", paste(CONTRAST_PRIMARY, collapse = " / "),
       " is present for ", TISSUE, ". Pick one from the table above and set ",
       "CONTRAST_PRIMARY.")
}
CONTRAST <- CONTRAST[1]
ok(paste0("selecting on contrast '", CONTRAST, "'"))

rna_c <- rna_t[contrast == CONTRAST]

subh("cell type tokens in the tidy table")
tok <- sort(unique(rna_c$celltype))
unmapped <- setdiff(tok, names(DISPLAY_LABELS))
for (tk in tok) {
  say(sprintf("  %-24s -> %s", tk,
              if (tk %in% names(DISPLAY_LABELS)) DISPLAY_LABELS[[tk]] else "NO MAPPING"))
}
if (length(unmapped)) {
  stop(length(unmapped), " cell type token(s) have no entry in DISPLAY_LABELS: ",
       paste(unmapped, collapse = ", "),
       "\nCopy them from DISPLAY_LABELS in multiome_common.py. Every figure ",
       "downstream joins on this, so a missing token is a silently wrong panel.")
}
ok(sprintf("all %d cell type tokens map to a display label", length(tok)))

# -----------------------------------------------------------------------------
# DUPLICATE CELL TYPE CHECK
#
# Added 2026-11-06 after the first brain full run selected VSIG4, NRP1, RPL29,
# PSMB6 and FOLR2 in BOTH Microglia-1 and Oligodendrocytes with byte-identical
# logFC and p-values. Two independent edgeR fits on different cells cannot do
# that, and microglia and oligodendrocytes have no biology in common that would
# explain it. One vendor file was read under two cell-type tokens.
#
# A duplicated cell type is worse than a missing one: it produces hits that
# look like independent replication of the same finding in a second population.
# So this runs before selection and refuses to continue.
# -----------------------------------------------------------------------------
bad_ct <- intersect(CORRUPTED_CELLTYPES[[TISSUE]], unique(rna_c$celltype))
if (length(bad_ct)) {
  subh("corrupted cell type labels, excluded")
  for (b in bad_ct) {
    say(sprintf("  %-24s dropped: its vendor table holds another cell type's numbers",
                b))
  }
  rna_c <- rna_c[!celltype %in% bad_ct]
  wrn(paste0(length(bad_ct), " cell type(s) excluded. No gene can be selected in ",
             "them and no figure should be drawn for them until Genohub ",
             "regenerates the tables."))
}

subh("duplicate cell type check")
wide <- dcast(rna_c[, .(celltype, gene, rna_logFC)], gene ~ celltype,
              value.var = "rna_logFC", fun.aggregate = function(x) x[1])
cts <- setdiff(names(wide), "gene")
dups <- data.table()
for (i in seq_along(cts)) {
  for (j in seq_len(i - 1L)) {
    va <- wide[[cts[i]]]; vb <- wide[[cts[j]]]
    k <- !is.na(va) & !is.na(vb)
    if (sum(k) < 50L) next
    if (isTRUE(all.equal(va[k], vb[k], tolerance = 0))) {
      dups <- rbind(dups, data.table(celltype_a = cts[i], celltype_b = cts[j],
                                     genes_compared = sum(k)))
    }
  }
}
if (nrow(dups) && !ALLOW_DUPLICATE_CELLTYPES) {
  print(dups)
  dump_csv(dups, sprintf("duplicate_celltypes_%s.csv", TISSUE))
  say("")
  stop(nrow(dups), " pair(s) of cell types carry IDENTICAL statistics across ",
       "every shared gene. That is one source file read under two tokens, not ",
       "biology.\n\n",
       "Fix it in MULTIOME_01_build_tidy.py before trusting any selection: a ",
       "duplicated cell type produces hits that look like independent ",
       "replication in a second population and are not.\n\n",
       "To look without fixing, set ALLOW_DUPLICATE_CELLTYPES <- TRUE, but do ",
       "not report anything from a run with that set.")
}
if (nrow(dups) && ALLOW_DUPLICATE_CELLTYPES) {
  print(dups)
  wrn("DUPLICATED CELL TYPES PRESENT and the override is on. This run is for")
  say("          diagnosis only. Do not report anything from it.")
} else {
  ok(sprintf("no duplicated cell types among %d compared pairwise", length(cts)))
}

# Cell types in the object with no RNA statistics at all. They can never be
# selected here, and ARC_05 may still be pooling them into a figure.
obj_ct <- OBJECT_CELLTYPES[[TISSUE]]
have_ct <- unname(DISPLAY_LABELS[tok])
no_rna <- setdiff(obj_ct, have_ct)
subh("cell types in the object with no RNA statistics")
if (length(no_rna)) {
  for (l in sort(no_rna)) say(paste0("  ", l))
  wrn(sprintf("%d of %d %s cell types have no RNA test, so no gene can be",
              length(no_rna), length(obj_ct), TISSUE))
  say("          selected in them. If ARC_05 pools a cell type by regex, check")
  say("          whether it is pulling in any of these.")
  if (TISSUE == "Brain" && any(grepl("^Microglia", no_rna))) {
    say("")
    wrn(paste0("ARC_05 pools ^Microglia, which includes ",
               paste(sort(grep("^Microglia", no_rna, value = TRUE)), collapse = " and "),
               "."))
    say("          Those cells appear in the coverage track but were never")
    say("          tested on the RNA side. Worth saying in the figure legend.")
  }
} else {
  ok("every cell type in the object has RNA statistics")
}

# -----------------------------------------------------------------------------
# Cell 6: Phase 0b, calibrate the sign of rna_logFC
#
# A contrast named VEH_vs_THCCBD does not tell us whether a positive logFC
# means higher in VEH or higher in THC/CBD; that depends on which group was
# ident.1 when the vendor ran it, and MULTIOME_01 left ATAC_IDENT1_GROUP unset
# for exactly this reason. The report gives us a calibration set: it states 22
# of these 23 genes as up in THC:CBD. If their logFC is predominantly positive,
# positive means up in THC:CBD. If predominantly negative, the contrast runs
# the other way and every sign in this script has to be flipped.
#
# This is evidence, not proof. It is reported with its agreement rate so you
# can see how strongly it leans before anything is read off a sign.
# -----------------------------------------------------------------------------

hdr("Phase 0b: sign convention of rna_logFC")

resolve_alias <- function(g, universe) {
  tries <- unique(c(g, GENE_ALIASES[[g]]))
  for (t in tries) if (t %in% universe) return(t)
  i <- which(toupper(universe) %in% toupper(tries))
  if (length(i)) return(universe[i[1]])
  NA_character_
}

gene_universe <- unique(rna_c$gene)
GENE_PANEL[, resolved := vapply(gene, resolve_alias, character(1), gene_universe)]

cal_genes <- GENE_PANEL[reported_rna %in% c("Up", "Slightly up") & !is.na(resolved), resolved]
cal <- rna_c[gene %in% cal_genes & !is.na(rna_logFC)]

SIGN_FLIP <- FALSE
if (nrow(cal) < 10) {
  wrn(sprintf("only %d rows available to calibrate the sign.", nrow(cal)))
  wrn("Direction is left unresolved rather than assumed.")
} else {
  # Weight by significance so near-zero noise does not dominate the vote.
  # rna_c is the source table, so this takes the source column name.
  cal_sig <- cal[get(if (FDR_BASIS == "glb") "rna_p_adj_glb" else "rna_p_adj_loc") < 0.05]
  use <- if (nrow(cal_sig) >= 10) cal_sig else cal
  pos <- mean(use$rna_logFC > 0)
  say(sprintf("  calibration set: %d gene x cell type rows (%s)",
              nrow(use), if (nrow(cal_sig) >= 10) "significant only" else "all"))
  say(sprintf("  the report calls these genes UP in THC:CBD"))
  say(sprintf("  fraction with positive rna_logFC: %.3f", pos))
  if (pos >= 0.65) {
    ok("positive rna_logFC means UP IN THC/CBD. No sign flip.")
  } else if (pos <= 0.35) {
    SIGN_FLIP <- TRUE
    wrn("positive rna_logFC means UP IN VEH. Signs are flipped for reporting.")
  } else {
    wrn(sprintf("the vote is not decisive (%.3f). Direction is NOT established.", pos))
    wrn("Signs are left as the table has them and every direction column is")
    wrn("marked unresolved. Ask the collaborators which group was ident.1.")
  }
}
SIGN_RESOLVED <- exists("pos") && (pos >= 0.65 || pos <= 0.35)
dir_label <- function(lfc) {
  if (!SIGN_RESOLVED) return(rep("unresolved", length(lfc)))
  up <- if (SIGN_FLIP) lfc < 0 else lfc > 0
  ifelse(is.na(lfc), NA_character_, ifelse(up, "up_in_THCCBD", "up_in_VEH"))
}

# -----------------------------------------------------------------------------
# Cell 7: Phase 1, provenance per gene and per cell type
# -----------------------------------------------------------------------------

hdr("Phase 1: provenance, per gene and per cell type")

celltypes <- sort(unique(rna_c$celltype))
grid <- CJ(gene = GENE_PANEL$gene, celltype = celltypes, unique = TRUE)
grid <- merge(grid, GENE_PANEL, by = "gene", all.x = TRUE)

# The raw p-value travels with every row. Without it an FDR of 0.88 cannot be
# told apart from "no effect" and "no power after correcting across 7,000
# genes", and those lead to completely different conversations.
rna_key <- rna_c[, .(celltype, gene,
                     rna_logFC,
                     rna_logCPM = if ("rna_logCPM" %in% names(rna_c)) rna_logCPM else NA_real_,
                     p_raw   = if ("rna_p_val" %in% names(rna_c)) rna_p_val else NA_real_,
                     fdr_loc = if ("rna_p_adj_loc" %in% names(rna_c)) rna_p_adj_loc else NA_real_,
                     fdr_glb = if ("rna_p_adj_glb" %in% names(rna_c)) rna_p_adj_glb else NA_real_,
                     has_linked_peak = if ("has_linked_peak" %in% names(rna_c)) has_linked_peak else NA)]
setnames(rna_key, "gene", "resolved")

prov <- merge(grid, rna_key, by = c("resolved", "celltype"), all.x = TRUE)

prov[, symbol_resolves := !is.na(resolved)]
prov[, tested_here := !is.na(rna_logFC)]
prov[, status := fifelse(!symbol_resolves, "symbol not in RNA table",
                  fifelse(!tested_here, "not detected or not tested in this cell type",
                          "tested"))]
prov[, celltype_display := unname(DISPLAY_LABELS[celltype])]

# Panel-local FDR: BH within each cell type across only the panel genes tested
# there. The 23 genes are a targeted list from the report, not a screen, so the
# genome-wide denominator behind fdr_loc answers a question nobody asked. This
# is a weaker claim and is labelled as such wherever it is used.
prov[, fdr_panel := NA_real_]
prov[tested_here == TRUE & !is.na(p_raw),
     fdr_panel := p.adjust(p_raw, method = "BH"), by = celltype]
prov[, n_panel_tested_in_celltype := sum(tested_here), by = celltype]

setcolorder(prov, c("gene", "resolved", "celltype", "celltype_display", "pathway",
                    "reported_rna", "reported_atac", "atac_possible", "weak_target",
                    "symbol_resolves", "tested_here", "status",
                    "rna_logFC", "p_raw", "fdr_loc", "fdr_glb", "fdr_panel"))

subh("symbol resolution against the RNA table")
sr <- unique(prov[, .(gene, resolved, symbol_resolves)])
for (i in seq_len(nrow(sr))) {
  r <- sr[i]
  if (!r$symbol_resolves) {
    wrn(sprintf("%-10s not found. Tried: %s", r$gene,
                paste(unique(c(r$gene, GENE_ALIASES[[r$gene]])), collapse = ", ")))
  } else if (r$resolved != r$gene) {
    ok(sprintf("%-10s found as '%s'", r$gene, r$resolved))
  }
}
ok(sprintf("%d of %d symbols resolve", sum(sr$symbol_resolves), nrow(sr)))

subh("coverage per gene across cell types")
cov <- prov[, .(tested_in = sum(tested_here),
                of = .N), by = .(gene, pathway)][order(pathway, -tested_in)]
print(cov)

dump_csv(prov, sprintf("gene_provenance_%s.csv", TISSUE))

# -----------------------------------------------------------------------------
# Cell 8: Phase 2, DEG selection within each cell type
# -----------------------------------------------------------------------------

hdr("Phase 2: RNA DEG selection, within each cell type")

hits <- prov[tested_here == TRUE]
hits[, abs_lfc := abs(rna_logFC)]
# Guard against the exact mistake that cost a run: a gate column name that is
# right for the source table and wrong for this one. Fails here, in one line,
# instead of inside a data.table eval frame.
if (!FDR_GATE_COL %in% names(hits)) {
  stop("FDR_GATE_COL '", FDR_GATE_COL, "' is not a column of hits. Present: ",
       paste(grep("^fdr|^p_", names(hits), value = TRUE), collapse = ", "))
}
hits[, fdr_gate := get(FDR_GATE_COL)]
hits[, passes := abs_lfc > LFC_PRIMARY & fdr_gate < FDR_PRIMARY]
hits[, direction := dir_label(rna_logFC)]

# What is actually binding? Printed before the sweep, because the first brain
# run failed entirely on FDR while every top effect cleared the logFC cut.
subh("what is binding")
n_lfc  <- hits[abs_lfc > LFC_PRIMARY, .N]
n_fdr  <- hits[fdr_gate < FDR_PRIMARY, .N]
say(sprintf("  rows tested:                       %d", nrow(hits)))
say(sprintf("  clearing |log2FC| > %-5.2f:          %d", LFC_PRIMARY, n_lfc))
say(sprintf("  clearing %s < %.2f:   %d", FDR_LABEL, FDR_PRIMARY, n_fdr))
say(sprintf("  clearing both:                     %d", hits[passes == TRUE, .N]))
say("")
say(sprintf("  best %s in the panel: %.4f", FDR_LABEL,
            suppressWarnings(min(hits$fdr_gate, na.rm = TRUE))))
say(sprintf("  best raw p in the panel:          %.3g",
            suppressWarnings(min(hits$p_raw, na.rm = TRUE))))
if (n_lfc > 0 && n_fdr == 0) {
  wrn("FDR is the binding constraint, not effect size. Lowering the logFC")
  say("          threshold cannot help. See the FDR comparison below.")
}

subh("the three corrections compared, at |log2FC| > 1")
big <- hits[abs_lfc > LFC_PRIMARY]
cmp <- data.table(
  correction = c("rna_p_adj_loc (within cluster, all genes)",
                 "rna_p_adj_glb (across all clusters)",
                 "fdr_panel (within cluster, panel genes only)"),
  best = c(suppressWarnings(min(big$fdr_loc, na.rm = TRUE)),
           suppressWarnings(min(big$fdr_glb, na.rm = TRUE)),
           suppressWarnings(min(big$fdr_panel, na.rm = TRUE))),
  n_under_05 = c(sum(big$fdr_loc < 0.05, na.rm = TRUE),
                 sum(big$fdr_glb < 0.05, na.rm = TRUE),
                 sum(big$fdr_panel < 0.05, na.rm = TRUE)),
  n_under_10 = c(sum(big$fdr_loc < 0.10, na.rm = TRUE),
                 sum(big$fdr_glb < 0.10, na.rm = TRUE),
                 sum(big$fdr_panel < 0.10, na.rm = TRUE)))
print(cmp)
say(sprintf("  selection is gated on FDR_BASIS = \"%s\"", FDR_BASIS))

# The panel is not independent of the data, so the panel-local correction has
# to be labelled every time it is used. Saying this once in a conversation is
# not enough; it has to travel with the output files.
if (FDR_BASIS == "panel") {
  say("")
  wrn("CIRCULARITY. The 23 genes come from the collaborators' report, which was")
  say("          derived from this same RNA dataset. A panel-local FDR computed")
  say("          on genes chosen because they already looked interesting in")
  say("          these data is not a valid error rate.")
  say("")
  say("          Legitimate: using it to PRIORITISE which ATAC tracks to draw.")
  say("          Not legitimate: reporting it as evidence that these genes are")
  say("          significantly differentially expressed.")
  say("")
  say("          The genome-wide number stays rna_p_adj_loc, carried in every")
  say("          output. The statistical weight belongs to the independent")
  say("          per-animal accessibility test in Phase 4.")
}

subh("threshold sweep, both axes")
sweep <- rbindlist(lapply(LFC_SWEEP, function(th) {
  rbindlist(lapply(FDR_SWEEP, function(fq) {
    x <- hits[abs_lfc > th & fdr_gate < fq]
    data.table(lfc_threshold = th, fdr_threshold = fq,
               gene_celltype_pairs = nrow(x),
               genes = uniqueN(x$gene),
               celltypes = uniqueN(x$celltype))
  }))
}))
print(dcast(sweep, lfc_threshold ~ fdr_threshold,
            value.var = "gene_celltype_pairs"))
say("  (rows are |log2FC| thresholds, columns are FDR thresholds, cells are")
say(sprintf("   gene x cell type pairs passing both, gated on %s)", FDR_LABEL))
dump_csv(sweep, sprintf("threshold_sweep_%s.csv", TISSUE))

sel <- hits[passes == TRUE][order(-abs_lfc)]
subh(sprintf("passing at |log2FC| > %.2f", LFC_PRIMARY))
if (!nrow(sel)) {
  wrn("NOTHING passes at the primary threshold.")
  say("")
  say("  This is the Gate B outcome. The sweep above says how much is there at")
  say("  0.25 and 0.50, and the per-gene best effect below says how far short")
  say("  each gene falls. That is what to take back to the collaborators.")
  say("")
  best <- hits[!is.na(abs_lfc), .SD[which.max(abs_lfc)],
               by = gene][order(-abs_lfc)]
  print(best[, .(gene, celltype_display, rna_logFC, p_raw, fdr_loc, fdr_panel,
                 direction)])
  say("")
  # What threshold would it actually take? Without this the obvious next move
  # is to lower the cut, and on the first brain run that could not have worked
  # at any defensible value.
  say("")
  say("  --- what threshold would be needed, on the current basis ---")
  bb <- hits[abs_lfc > LFC_PRIMARY & !is.na(fdr_gate)][order(fdr_gate)]
  if (nrow(bb)) {
    for (k in c(1L, 3L, 5L, 10L)) {
      if (k <= nrow(bb)) {
        say(sprintf("    to capture %2d pair(s): FDR < %.3f   (%s in %s)", k,
                    bb$fdr_gate[k], bb$gene[k], bb$celltype_display[k]))
      }
    }
    if (bb$fdr_gate[1] > 0.25) {
      say("")
      wrn("Even one pair needs an FDR above 0.25. Lowering the cut is not the")
      say("          answer here; the denominator is. See the correction table")
      say("          above and consider FDR_BASIS = \"panel\".")
    }
  }
  say("")
  say("  Read p_raw against fdr_loc. A small p_raw with a large fdr_loc is a")
  say("  multiple-testing result and the panel-local column is the fairer test.")
  say("  A large p_raw is no effect, and no correction rescues it.")
  dump_csv(best, sprintf("best_effect_per_gene_%s.csv", TISSUE))
} else {
  print(sel[, .(gene, celltype_display, rna_logFC, p_raw, fdr_gate, direction,
                pathway)])
  ok(sprintf("%d gene x cell type pairs, %d genes, %d cell types",
             nrow(sel), uniqueN(sel$gene), uniqueN(sel$celltype)))
  if (any(sel$weak_target)) {
    wrn(sprintf("%d of these are housekeeping loci (%s) whose promoters are open",
                sum(sel$weak_target),
                paste(unique(sel[weak_target == TRUE, gene]), collapse = ", ")))
    say("          in nearly every cell type. An accessibility result there is")
    say("          close to guaranteed and says little.")
  }
  if (any(!sel$atac_possible)) {
    wrn(sprintf("%d pass on RNA but can never have an ATAC panel: %s",
                sum(!sel$atac_possible),
                paste(unique(sel[atac_possible == FALSE, gene]), collapse = ", ")))
  }
}
# Where the data disagrees with the report's stated direction. The report
# calls 22 of 23 genes up in THC:CBD; a gene whose best effect runs the other
# way is a discrepancy to raise, not a result to quietly plot.
if (SIGN_RESOLVED) {
  subh("direction against the report")
  bestd <- hits[!is.na(abs_lfc), .SD[which.max(abs_lfc)], by = gene]
  bestd <- merge(bestd, GENE_PANEL[, .(gene, reported_rna)], by = "gene",
                 suffixes = c("", ".p"))
  bestd[, reported_up := reported_rna %in% c("Up", "Slightly up")]
  bestd[, observed_up := direction == "up_in_THCCBD"]
  disagree <- bestd[reported_up != observed_up]
  if (nrow(disagree)) {
    wrn(sprintf("%d gene(s) run opposite to the report at their strongest cell type:",
                nrow(disagree)))
    print(disagree[, .(gene, reported_rna, celltype_display, rna_logFC,
                       p_raw, direction)])
  } else {
    ok("every gene agrees with the report's stated direction")
  }
}

# -----------------------------------------------------------------------------
# DROP LEDGER
#
# Every gene x cell type pair gets exactly one reason. Without this, "the gene
# is not in the output" covers five completely different situations that call
# for five different responses, and the collaborators cannot tell them apart.
# -----------------------------------------------------------------------------
ledger <- copy(prov)
ledger[, drop_reason := fifelse(
  !symbol_resolves,                        "absent from the RNA table entirely",
  fifelse(!tested_here,                    "not expressed or not tested in this cell type",
  fifelse(abs(rna_logFC) <= LFC_PRIMARY,   sprintf("effect below |log2FC| %.2f", LFC_PRIMARY),
  fifelse(get(FDR_GATE_COL) >= FDR_PRIMARY, sprintf("FDR above %.2f on the %s basis", FDR_PRIMARY, FDR_BASIS),
                                           "selected on RNA"))))]
ledger[symbol_resolves & !atac_possible,
       drop_reason := paste0(drop_reason, "; and can never have an ATAC panel (",
                             atac_impossible_reason, ")")]
dump_csv(ledger, sprintf("drop_ledger_%s.csv", TISSUE))

subh("why each gene is or is not in the shortlist")
led <- ledger[, .(cell_types = .N), by = .(gene, drop_reason)]
# One row per gene: its best showing anywhere. A gene selected in any cell type
# reports as selected; otherwise its strongest effect explains why not.
ledger[, rank_key := fifelse(drop_reason == "selected on RNA", 0L, 1L)]
best_reason <- ledger[order(rank_key, -abs(rna_logFC), na.last = TRUE),
                      .SD[1], by = gene]
print(best_reason[order(rank_key, -abs(rna_logFC), na.last = TRUE),
                  .(gene, pathway, celltype_display, rna_logFC, p_raw,
                    drop_reason)])

subh("drop reasons, counted by gene")
print(ledger[, .(gene_celltype_pairs = .N,
                 genes = uniqueN(gene)), by = drop_reason][order(-gene_celltype_pairs)])

dump_csv(hits, sprintf("deg_hits_%s.csv", TISSUE))

# --- the matrix: cell types down, genes across ------------------------------

subh("cell type by gene matrix")
mat <- dcast(hits, celltype_display ~ gene, value.var = "rna_logFC")
dump_csv(mat, sprintf("deg_matrix_%s.csv", TISSUE))

mat_fdr <- dcast(hits, celltype_display ~ gene, value.var = "fdr_loc")
dump_csv(mat_fdr, sprintf("deg_matrix_fdr_%s.csv", TISSUE))

hm <- copy(hits)
hm[, star := fifelse(passes, "*", "")]
gene_order <- GENE_PANEL[order(pathway, gene), gene]
hm[, gene := factor(gene, levels = intersect(gene_order, unique(hm$gene)))]
ct_order <- hits[, .(m = max(abs_lfc, na.rm = TRUE)), by = celltype_display][order(-m), celltype_display]
hm[, celltype_display := factor(celltype_display, levels = rev(ct_order))]

lim <- suppressWarnings(max(abs(hm$rna_logFC), na.rm = TRUE))
if (!is.finite(lim) || lim <= 0) lim <- 1
p <- ggplot(hm, aes(x = gene, y = celltype_display, fill = rna_logFC)) +
  geom_tile(color = "#FFFFFF", linewidth = 0.4) +
  # geom_text size is in mm, not points: pt / 2.845 keeps the star at the
  # same visual weight as FONT_SIZE elsewhere.
  geom_text(aes(label = star), size = FONT_SIZE / 2.845, vjust = 0.78) +
  scale_fill_gradient2(low = HEAT_LOW, mid = HEAT_MID, high = HEAT_HIGH,
                       midpoint = 0, limits = c(-lim, lim), na.value = "#E8E8E8",
                       name = "log2FC") +
  labs(title = sprintf("%s: report gene panel by cell type", TISSUE),
       subtitle = sprintf("contrast %s  |  * marks |log2FC| > %.2f and %s < %.2f  |  grey = not tested in that cell type",
                          CONTRAST, LFC_PRIMARY, FDR_LABEL, FDR_PRIMARY),
       x = NULL, y = NULL) +
  theme_minimal(base_size = FONT_SIZE) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = FONT_SIZE * 0.72, color = "#5A5A5A"),
        panel.grid = element_blank())

fig_w <- 8 + 0.85 * uniqueN(hm$gene)
fig_h <- 5 + 0.70 * uniqueN(hm$celltype_display)
stem <- file.path(OUT_DIR, sprintf("deg_matrix_%s", TISSUE))
ggsave(paste0(stem, ".pdf"), p, width = fig_w, height = fig_h,
       units = "in", device = grDevices::cairo_pdf, limitsize = FALSE)
ggsave(paste0(stem, ".png"), p, width = fig_w, height = fig_h,
       units = "in", dpi = FIG_DPI, limitsize = FALSE)
say(sprintf("    wrote deg_matrix_%s.pdf and .png  (%.0f x %.0f in, %d dpi)",
            TISSUE, fig_w, fig_h, FIG_DPI))

if (MODE == "select") {
  hdr("Stopping after Phase 2 (select mode)")
  say(sprintf("  Read deg_matrix_%s.png first, then:", TISSUE))
  say(sprintf("    Rscript ARC_06_select_celltypes.R %s full", TISSUE))
  say("")
  say(sprintf("  elapsed: %.1f min",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  quit(save = "no", status = 0)
}

if (!nrow(sel)) {
  hdr("Stopping: nothing passed Gate B")
  say("  There is no shortlist to inventory peaks for. Take the sweep and the")
  say("  per-gene best effects above back to the collaborators first.")
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Cell 9: Phase 3, peak inventory within the surviving cell types
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(GenomicRanges); library(Matrix)
})

hdr("Phase 3: peak inventory, within the surviving cell types")

rds <- file.path(OBJ_ROOT, TISSUE, sprintf("seurat_%s_rebuilt.rds", TISSUE))
if (!file.exists(rds)) stop("not found: ", rds, "\nRun ARC_04_build_object.R first.")
say(paste0("  reading ", rds))
obj <- readRDS(rds)
DefaultAssay(obj) <- "ATAC"
ok(sprintf("%s cells, %s peaks",
           format(ncol(obj), big.mark = ","),
           format(nrow(obj[["ATAC"]]), big.mark = ",")))

# Gate on the join before using it. This is the whole reason DISPLAY_LABELS is
# allowed to be duplicated in this file.
obj_labels <- unique(as.character(obj$celltype_display))
want_labels <- unique(sel$celltype_display)
missing_lab <- setdiff(want_labels, obj_labels)
if (length(missing_lab)) {
  say("")
  say("  labels in the object:")
  for (l in sort(obj_labels)) say(paste0("    ", l))
  stop(length(missing_lab), " selected cell type(s) are not in the object: ",
       paste(missing_lab, collapse = ", "),
       "\nDISPLAY_LABELS in this script disagrees with the object. Fix it here.")
}
ok(sprintf("all %d selected cell type labels found in the object", length(want_labels)))

ann <- Annotation(obj[["ATAC"]])
pk <- StringToGRanges(rownames(obj[["ATAC"]]), sep = c("-", "-"))
counts <- GetAssayData(obj, assay = "ATAC", layer = "counts")

obj$animal <- as.character(obj$animal)
obj$arm <- unname(ANIMAL_ARM[obj$animal])
# Post-infection only, matching the arm comparison ARC_05 draws.
post <- if ("timepoint" %in% names(obj@meta.data)) {
  as.character(obj$timepoint) != "Pre"
} else rep(TRUE, ncol(obj))

gene_window <- function(symbol) {
  g <- ann[!is.na(ann$gene_name) & ann$gene_name == symbol]
  if (!length(g)) return(NULL)
  gr <- GRanges(as.character(seqnames(g))[1],
                IRanges(min(start(g)), max(end(g))),
                strand = as.character(strand(g))[1])
  tss <- if (as.character(strand(gr)) == "-") end(gr) else start(gr)
  w <- gr
  if (as.character(strand(gr)) == "-") {
    start(w) <- max(1, start(gr) - WINDOW_DOWN); end(w) <- end(gr) + WINDOW_UP
  } else {
    start(w) <- max(1, start(gr) - WINDOW_UP);   end(w) <- end(gr) + WINDOW_DOWN
  }
  list(gene = gr, window = w, tss = tss)
}

inv <- list()
for (i in seq_len(nrow(sel))) {
  r <- sel[i]
  if (!isTRUE(r$atac_possible)) next
  gw <- gene_window(r$resolved)
  if (is.null(gw)) {
    wrn(sprintf("%s: not in the ATAC annotation, skipped", r$gene)); next
  }
  idx <- which(countOverlaps(pk, gw$window) > 0)
  if (!length(idx)) {
    wrn(sprintf("%s / %s: no consensus peaks in the window", r$gene, r$celltype_display))
    next
  }
  in_ct <- post & as.character(obj$celltype_display) == r$celltype_display
  if (!any(in_ct)) next

  mid <- (start(pk[idx]) + end(pk[idx])) %/% 2
  dist <- mid - gw$tss
  if (as.character(strand(gw$gene)) == "-") dist <- -dist

  sub_counts <- counts[idx, in_ct, drop = FALSE]
  arm_v <- obj$arm[in_ct]
  # vapply, not sapply: with a single peak in the window sapply returns a
  # bare vector and per_arm[, "VEH"] then fails.
  per_arm <- vapply(c("VEH", "THC/CBD"), function(a) {
    k <- which(arm_v == a)
    if (!length(k)) return(rep(NA_real_, length(idx)))
    as.numeric(Matrix::rowMeans(sub_counts[, k, drop = FALSE]))
  }, numeric(length(idx)))
  if (is.null(dim(per_arm))) per_arm <- matrix(per_arm, nrow = 1,
                                               dimnames = list(NULL, c("VEH", "THC/CBD")))

  d <- data.table(
    gene = r$gene, resolved = r$resolved,
    celltype = r$celltype, celltype_display = r$celltype_display,
    peak = rownames(obj[["ATAC"]])[idx],
    chrom = as.character(seqnames(pk[idx])),
    start = start(pk[idx]), end = end(pk[idx]), width = width(pk[idx]),
    dist_to_tss = dist,
    region_class = fifelse(abs(dist) <= PROMOTER_BP, "promoter",
                    fifelse(dist < 0, "upstream_distal", "downstream_distal")),
    cells_VEH = sum(arm_v == "VEH"), cells_THCCBD = sum(arm_v == "THC/CBD"),
    mean_count_VEH = per_arm[, "VEH"], mean_count_THCCBD = per_arm[, "THC/CBD"]
  )
  if (!is.null(da)) {
    dag <- da[tissue == TISSUE & celltype == r$celltype]
    if (nrow(dag)) {
      dgr <- GRanges(dag$chrom, IRanges(dag$start, dag$end))
      d[, vendor_da_peak := countOverlaps(pk[idx], dgr) > 0]
    }
  }
  inv[[length(inv) + 1]] <- d
  say(sprintf("  %-10s %-28s %3d peaks (%d promoter)  VEH n=%d  THC/CBD n=%d",
              r$gene, r$celltype_display, nrow(d),
              sum(d$region_class == "promoter"), d$cells_VEH[1], d$cells_THCCBD[1]))
}

peak_support <- if (length(inv)) rbindlist(inv, fill = TRUE) else data.table()
if (!nrow(peak_support)) {
  hdr("Stopping: no peaks to test")
  quit(save = "no", status = 0)
}
dump_csv(peak_support, sprintf("peak_support_%s.csv", TISSUE))

# -----------------------------------------------------------------------------
# Cell 10: Gate C and Phase 4, per-animal pseudobulk within cell type
#
# edgeR quasi-likelihood on counts summed per animal within the cell type. This
# is the same model muscat's pbDS() wraps; it is run directly here because the
# input is a peak subset rather than a full SingleCellExperiment, and going
# through muscat would mean converting the whole object to test a few hundred
# features.
# -----------------------------------------------------------------------------

hdr("Phase 4: per-animal pseudobulk DA, within cell type")

have_edger <- requireNamespace("edgeR", quietly = TRUE)
if (!have_edger) {
  wrn("edgeR is not installed, so Phase 4 is skipped. Install bioconductor-edger.")
}

da_rows <- list()
for (ctd in unique(peak_support$celltype_display)) {
  in_ct <- post & as.character(obj$celltype_display) == ctd
  an <- obj$animal[in_ct]
  arm_v <- obj$arm[in_ct]
  per_an <- table(an)
  keep_an <- names(per_an)[per_an >= MIN_CELLS_PER_ANIMAL]
  arm_of <- ANIMAL_ARM[keep_an]
  n_veh <- sum(arm_of == "VEH"); n_thc <- sum(arm_of == "THC/CBD")

  subh(sprintf("%s", ctd))
  for (a in names(sort(per_an, decreasing = TRUE))) {
    say(sprintf("  %-8s %-9s %5d cells%s", a, ANIMAL_ARM[[a]], per_an[[a]],
                if (per_an[[a]] < MIN_CELLS_PER_ANIMAL) "   (below floor)" else ""))
  }
  say(sprintf("  animals clearing the %d-cell floor: VEH %d, THC/CBD %d",
              MIN_CELLS_PER_ANIMAL, n_veh, n_thc))

  if (!have_edger || n_veh < MIN_ANIMALS_PER_ARM || n_thc < MIN_ANIMALS_PER_ARM) {
    wrn(sprintf("Gate C not met (need %d per arm). No statistic for this cell type.",
                MIN_ANIMALS_PER_ARM))
    say("          The coverage figure is still drawn, labelled descriptive only.")
    next
  }

  idx <- match(unique(peak_support[celltype_display == ctd, peak]),
               rownames(obj[["ATAC"]]))
  idx <- idx[!is.na(idx)]
  cols <- which(in_ct)
  pb <- vapply(keep_an, function(a) {
    k <- cols[an == a]
    as.numeric(Matrix::rowSums(counts[idx, k, drop = FALSE]))
  }, numeric(length(idx)))
  if (is.null(dim(pb))) pb <- matrix(pb, nrow = length(idx),
                                     dimnames = list(NULL, keep_an))
  rownames(pb) <- rownames(obj[["ATAC"]])[idx]

  grp <- factor(unname(ANIMAL_ARM[colnames(pb)]), levels = c("VEH", "THC/CBD"))
  y <- edgeR::DGEList(counts = pb, group = grp)
  y <- y[edgeR::filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
  if (!nrow(y)) { wrn("no peaks survive filterByExpr here"); next }
  y <- edgeR::calcNormFactors(y)
  des <- model.matrix(~ grp)
  y <- edgeR::estimateDisp(y, des)
  fit <- edgeR::glmQLFit(y, des)
  res <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2), n = Inf)$table

  d <- as.data.table(res, keep.rownames = "peak")
  d[, `:=`(celltype_display = ctd,
           n_animals_VEH = n_veh, n_animals_THCCBD = n_thc,
           model = "edgeR QL, per-animal pseudobulk, within cell type")]
  da_rows[[length(da_rows) + 1]] <- d
  ok(sprintf("%d peaks tested, %d at FDR < 0.05", nrow(d), sum(d$FDR < 0.05)))
}

da_pb <- if (length(da_rows)) rbindlist(da_rows, fill = TRUE) else data.table()
if (nrow(da_pb)) dump_csv(da_pb, sprintf("da_pseudobulk_%s.csv", TISSUE))

# -----------------------------------------------------------------------------
# Cell 11: Phase 5, emit the ARC_05 panel
# -----------------------------------------------------------------------------

hdr("Phase 5: ARC_05 panel")

panel <- unique(peak_support[, .(gene, resolved, celltype, celltype_display)])
panel[, tissue := TISSUE]
# ARC_05 matches cell type by regex against celltype_display. Anchored and
# escaped so "Microglia-1" cannot also match "Microglia-10" if one ever exists.
panel[, celltype_regex := paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1",
                                           celltype_display), "$")]
panel[, celltype_label := celltype_display]
panel[, window_up := WINDOW_UP]
panel[, window_down := WINDOW_DOWN]
panel <- merge(panel, GENE_PANEL[, .(gene, pathway, weak_target)], by = "gene", all.x = TRUE)
if (nrow(da_pb)) {
  best <- da_pb[, .(min_da_fdr = min(FDR, na.rm = TRUE)), by = celltype_display]
  panel <- merge(panel, best, by = "celltype_display", all.x = TRUE)
  panel[, statistic := fifelse(is.na(min_da_fdr), "descriptive only", "pseudobulk DA")]
} else {
  panel[, statistic := "descriptive only"]
}
panel[, selection_basis := FDR_BASIS]
panel[, selection_note := fifelse(
  FDR_BASIS == "panel",
  "prioritised on a panel-local FDR; panel not independent of the data; not a reported error rate",
  "selected on the genome-wide within-cluster FDR")]
setcolorder(panel, c("tissue", "gene", "resolved", "celltype_label",
                     "celltype_regex", "window_up", "window_down",
                     "pathway", "weak_target", "statistic",
                     "selection_basis", "selection_note"))
print(panel)
dump_csv(panel, sprintf("panel_%s.csv", TISSUE))

hdr("Done")
say(sprintf("  %d gene x cell type pair(s) selected", nrow(panel)))
say(sprintf("  output: %s", OUT_DIR))
say(sprintf("  elapsed: %.1f min",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
say("")
say("  Next: point ARC_05 at panel_<tissue>.csv.")
if (FDR_BASIS == "panel") {
  say("")
  wrn("These pairs were PRIORITISED, not shown to be significant. Any figure")
  say("          legend built from them should cite the Phase 4 pseudobulk FDR,")
  say("          not the RNA selection.")
}
if (!SIGN_RESOLVED) {
  say("")
  wrn("The sign of rna_logFC was NOT established in Phase 0b, so every")
  say("          direction column reads 'unresolved'. Do not state up or down")
  say("          anywhere until the collaborators say which group was ident.1.")
}
invisible(gc(FALSE))
