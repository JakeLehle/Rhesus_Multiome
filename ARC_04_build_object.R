#!/usr/bin/env Rscript
# =============================================================================
# ARC_04_build_object.R  <tissue>
#
# Builds one Seurat object per tissue from the cellranger-arc count output, with
# the vendor's cell annotations transferred, and a ChromatinAssay carrying the
# per-library fragment files so Signac can draw coverage tracks.
#
#   Rscript ARC_04_build_object.R Brain
#   Rscript ARC_04_build_object.R Colon
#
# Normally submitted through ARC_04_submit.sh.
#
# -----------------------------------------------------------------------------
# DECISIONS THIS SCRIPT ENCODES, AND WHERE THEY CAME FROM
# -----------------------------------------------------------------------------
#
# CELL NAMES ARE SET EXPLICITLY, NOT BY merge().
#   The vendor's CellID is BARCODE-1_<index>: a 16 bp barcode, the -1 that every
#   cellranger filtered matrix carries, and the list position Seurat's merge()
#   appends. ARC_00 recovered the index order and it is the alphabetical sort of
#   the loom Samples labels.
#
#   We do NOT rely on merge() to produce those suffixes. merge() only disambiguates
#   names that actually collide, and only 137 brain and 537 colon barcodes are
#   reused across libraries, so it would suffix some cells and leave others bare.
#   Setting the names ourselves makes every cell match the vendor's exactly.
#
# THE ANNOTATION JOIN IS AN INNER JOIN.
#   ARC_03 measured retention: brain median 75.5%, colon median 29.7%. The vendor
#   applied their own QC after count, roughly nUMI >= 1,000 and nGene >= 300. An
#   outer join would put tens of thousands of unannotated cells into the coverage
#   plots, and they would not be the cells behind the collaborators' figures.
#
# THE CONSENSUS PEAK SET IS A REDUCED UNION OF THE PER-LIBRARY BEDs.
#   ARC_03 settled this. The union gives 203,193 peaks for brain against the
#   vendor's 204,136 (ratio 0.995) and 339,504 for colon against 345,553 (0.982),
#   and 99.9% of their DA peak midpoints fall inside it. So `cellranger-arc aggr`
#   is unnecessary, which also means no barcode-suffix remapping.
#
# FEATURES ARE cbind-ed, NOT MERGED.
#   Every library was counted against the same Mmul8 reference, so the GEX
#   feature set is identical across all of them. A direct sparse cbind is faster
#   and less memory-hungry than Seurat merge(), and the script verifies the
#   rownames match before doing it.
#
# FeatureMatrix IS CHECKPOINTED PER LIBRARY.
#   It is the expensive step: ~340K peaks against 8 or 12 fragment files. Each
#   library's result is written to disk, so a failure or a timeout resumes rather
#   than restarting.
#
# Jake Lehle, 2026-10
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Cell 1: parameters
# -----------------------------------------------------------------------------

ANALYSIS_ROOT <- "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
COUNT_ROOT    <- file.path(ANALYSIS_ROOT, "02_arc", "02_count")
LOOM_META_DIR <- file.path(ANALYSIS_ROOT, "02_arc", "00_loom_metadata")
OUT_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "04_objects")

GTF <- "/master/jlehle/WORKING/REFS/Mmul8/genes.gtf"

# Library order per tissue: the vendor's Seurat merge() position, recovered by
# ARC_00. Position in this vector IS the _<index> suffix.
LIB_ORDER <- list(
  Brain = c("40872", "40873", "40874", "40875", "40876", "40877", "40878", "40879"),
  Colon = c("40872N", "40873N", "40874N", "40874P", "40875N", "40875P",
            "40876N", "40877N", "40877P", "40878N", "40878P", "40879N")
)

# Consensus peak filtering. Standard contigs only; chrM is excluded because the
# mkref config declared it non-nuclear and the vendor's data has no chrM peaks.
KEEP_CONTIGS <- c(paste0("chr", 1:20), "chrX", "chrY")
PEAK_MIN_WIDTH <- 50
PEAK_MAX_WIDTH <- 10000

# FeatureMatrix chunking. Lower if memory becomes a problem.
FM_CHUNK <- 10000

# CreateFragmentObject validates that the requested cells appear in the file,
# which reads through it. Slow but worth it once; set FALSE on a rerun.
VALIDATE_FRAGMENTS <- TRUE

# Reuse the per-library FeatureMatrix checkpoints if present.
USE_CHECKPOINTS <- TRUE

# -----------------------------------------------------------------------------
# Cell 2: plumbing
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: ARC_04_build_object.R <Brain|Colon>")
TISSUE <- args[1]
if (!TISSUE %in% names(LIB_ORDER)) stop("tissue must be Brain or Colon")

OUT_DIR <- file.path(OUT_ROOT, TISSUE)
CKPT_DIR <- file.path(OUT_DIR, "featurematrix_checkpoints")
dir.create(CKPT_DIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(paste0(..., "\n"))
hdr <- function(x) { say(""); say(strrep("=", 76)); say(x); say(strrep("=", 76)) }
# NOT named sub(): that shadows base::sub, which this script needs for the
# sample-suffix and cell-name regexes further down.
subh <- function(x) { say(""); say(paste0("--- ", x, " ---")) }
ok  <- function(x) say(paste0("  [ok]    ", x))
wrn <- function(x) say(paste0("  [warn]  ", x))

# Regression guard for exactly that mistake. It checks this script's own helper
# names against base directly, rather than testing whether sub() still resolves
# to base::sub, because Bioconductor turns several base functions into S4
# generics and that test would false-alarm. A shadowed base function fails 20
# minutes in, after the GEX read, so it is worth catching at startup.
for (fn in c("say", "hdr", "subh", "ok", "wrn", "read_gex_matrix",
             "call_supported")) {
  if (exists(fn, envir = baseenv(), inherits = FALSE)) {
    stop("helper '", fn, "' collides with base::", fn,
         ". Rename it before running, or calls to the base function will ",
         "silently dispatch to the helper.")
  }
}

# Signac's signatures move between releases. validate.cells was dropped from
# CreateFragmentObject by 1.17, and because that function takes ..., an
# unsupported argument is not rejected at the call: it falls through into
# ValidateCells and errors there as "unused argument", 15 minutes into a run and
# after the fragment hash has already been computed.
#
# So optional arguments are filtered against the function's own formals first.
# An argument that this version does not take is dropped, with a line saying so,
# and the function's default applies. Required arguments are passed through
# untouched, because a missing required argument should still be a hard error.
call_supported <- function(fn, required, optional = list(), label = "call") {
  f <- names(formals(fn))
  keep <- optional[names(optional) %in% f]
  dropped <- setdiff(names(optional), names(keep))
  if (length(dropped)) {
    wrn(paste0(label, ": this version does not take ",
               paste(dropped, collapse = ", "),
               "; the function's own default applies"))
  }
  do.call(fn, c(required, keep))
}

t_start <- Sys.time()
hdr(paste0("ARC_04_build_object : ", TISSUE))
say(paste0("host:      ", Sys.info()[["nodename"]]))
say(paste0("started:   ", format(t_start, "%Y-%m-%d %H:%M:%S")))
say(paste0("out:       ", OUT_DIR))
say(paste0("libraries: ", length(LIB_ORDER[[TISSUE]])))
say("")
say("package versions:")
for (p in c("Seurat", "Signac", "Matrix", "GenomicRanges")) {
  say(paste0("  ", formatC(p, width = -16), as.character(packageVersion(p))))
}

libs <- LIB_ORDER[[TISSUE]]
run_ids <- paste0(TISSUE, "_", libs)
names(run_ids) <- seq_along(run_ids)   # index -> run_id

# -----------------------------------------------------------------------------
# Cell 3: the target cell list, from the loom metadata
# -----------------------------------------------------------------------------

hdr("Phase 1: target cells from the loom")

meta_path <- file.path(LOOM_META_DIR, paste0("cell_metadata_", TISSUE, ".csv"))
if (!file.exists(meta_path)) stop("not found: ", meta_path, "\nRun ARC_00_loom_metadata.py first.")

md <- fread(meta_path)
say(paste0("  loom cells: ", format(nrow(md), big.mark = ",")))

req <- c("CellID", "barcode", "library_index", "Samples", "Stimulation",
         "celltype", "celltype_display", "umap_1", "umap_2")
missing_cols <- setdiff(req, names(md))
if (length(missing_cols)) {
  stop("cell_metadata is missing columns: ", paste(missing_cols, collapse = ", "))
}

# The vendor's loom carries literal newlines inside some cluster labels, e.g.
# "Entero\ncyte-3", which is why a wc -l of this CSV can exceed its row count.
# A label with a newline in it survives into the object and then into every
# figure legend and filename built from it, so they are collapsed here and the
# count is reported rather than fixed silently.
label_cols <- c("Samples", "Stimulation", "celltype", "celltype_display")
n_fixed <- 0L
for (cl in label_cols) {
  v <- as.character(md[[cl]])
  w2 <- gsub("[\r\n]+", " ", v)
  w2 <- trimws(gsub("[[:space:]]{2,}", " ", w2))
  n_fixed <- n_fixed + sum(w2 != v, na.rm = TRUE)
  md[[cl]] <- w2
}
if (n_fixed > 0) {
  wrn(paste0(n_fixed, " label value(s) contained embedded newlines or padded ",
             "whitespace and were collapsed. Check ARC_00's clean_label() if ",
             "you expected it to have handled this upstream."))
} else {
  ok("no embedded newlines in the label columns")
}

md$library_index <- as.integer(md$library_index)
stopifnot(all(md$library_index %in% seq_along(libs)))

subh("cells per library")
tab <- md[, .N, by = .(library_index, Samples)][order(library_index)]
for (i in seq_len(nrow(tab))) {
  say(sprintf("  _%-3d %-10s %7s cells",
              tab$library_index[i], tab$Samples[i],
              format(tab$N[i], big.mark = ",")))
}

# The index-to-sample mapping must agree with LIB_ORDER or every cell name is
# wrong and the join silently produces nothing.
for (i in seq_along(libs)) {
  s <- unique(md[library_index == i, Samples])
  if (length(s) != 1 || s != libs[i]) {
    stop("library index ", i, " maps to sample '", paste(s, collapse = ","),
         "' but LIB_ORDER expects '", libs[i], "'. ",
         "The ordering recovered by ARC_00 and the one hardcoded here disagree; ",
         "fix LIB_ORDER before going further.")
  }
}
ok("library index to sample mapping agrees with LIB_ORDER")

subh("cell types")
ct <- md[, .N, by = celltype_display][order(-N)]
for (i in seq_len(nrow(ct))) {
  say(sprintf("  %-32s %7s", ct$celltype_display[i], format(ct$N[i], big.mark = ",")))
}

# -----------------------------------------------------------------------------
# Cell 4: GEX matrices
#
# One cbind rather than a Seurat merge. Every library used the same reference so
# the features are identical, which the check below enforces.
# -----------------------------------------------------------------------------

hdr("Phase 2: gene expression")

# cellranger-arc writes the same filtered matrix twice: as an HDF5 file and as
# an MEX directory. We prefer the directory, because Read10X goes through
# Matrix::readMM and needs nothing beyond base R, while Read10X_h5 needs hdf5r.
# That is one less package to install, and the two routes return the same
# matrix. strip.suffix = FALSE is explicit: the trailing -1 is part of the
# vendor's CellID and dropping it would break every join downstream.
read_gex_matrix <- function(rid) {
  mex <- file.path(COUNT_ROOT, rid, "outs", "filtered_feature_bc_matrix")
  h5  <- file.path(COUNT_ROOT, rid, "outs", "filtered_feature_bc_matrix.h5")

  mex_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  if (dir.exists(mex) && all(file.exists(file.path(mex, mex_files)))) {
    return(list(x = Read10X(mex, strip.suffix = FALSE), via = "MEX"))
  }

  if (file.exists(h5)) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
      stop(rid, ": the MEX directory is incomplete and reading the .h5 instead ",
           "needs hdf5r, which is not installed. Either restore\n  ", mex,
           "\nor install r-hdf5r into this environment.")
    }
    return(list(x = Read10X_h5(h5), via = "HDF5"))
  }

  stop(rid, ": no filtered matrix found. Looked for\n  ", mex, "\n  ", h5)
}

gex_list <- vector("list", length(libs))
feat_ref <- NULL
match_report <- data.table()
read_via <- character(0)

for (i in seq_along(libs)) {
  rid <- run_ids[i]

  r10 <- read_gex_matrix(rid)
  x <- r10$x
  read_via <- c(read_via, r10$via)
  gex <- if (is.list(x)) x[["Gene Expression"]] else x
  if (is.null(gex)) stop(rid, ": no 'Gene Expression' matrix in the ", r10$via, " input")

  # The explicit rename. Barcodes already end in -1, so appending _<index>
  # reproduces the vendor's CellID exactly.
  colnames(gex) <- paste0(colnames(gex), "_", i)

  if (is.null(feat_ref)) {
    feat_ref <- rownames(gex)
  } else if (!identical(rownames(gex), feat_ref)) {
    stop(rid, ": feature set differs from the first library. Every library must ",
         "be counted against the same reference or a cbind is invalid.")
  }

  want <- md[library_index == i, CellID]
  keep <- intersect(want, colnames(gex))
  gex_list[[i]] <- gex[, keep, drop = FALSE]

  match_report <- rbind(match_report, data.table(
    index = i, run_id = rid, sample = libs[i],
    cells_in_count = ncol(gex), cells_in_loom = length(want),
    matched = length(keep),
    pct_loom_found = round(100 * length(keep) / max(length(want), 1), 2)
  ))

  say(sprintf("  _%-3d %-14s %-4s  count %6d   loom %6d   matched %6d  (%5.1f%% of loom)",
              i, rid, r10$via, ncol(gex), length(want), length(keep),
              100 * length(keep) / max(length(want), 1)))
  rm(gex, x, r10); invisible(gc(FALSE))
}

if (length(unique(read_via)) > 1) {
  wrn(paste0("matrices were read by more than one route (",
             paste(unique(read_via), collapse = ", "),
             "), which is fine but means some MEX directories are incomplete"))
}

subh("barcode match rate")
print(match_report)

# The acceptance gate. Anything short of a near-complete match means the cell
# naming convention is wrong, and every downstream join would be quietly empty.
bad <- match_report[pct_loom_found < 99]
if (nrow(bad)) {
  say("")
  for (i in seq_len(nrow(bad))) {
    wrn(sprintf("%s: only %.2f%% of the loom's cells were found in the count output",
                bad$run_id[i], bad$pct_loom_found[i]))
  }
  stop("barcode matching failed for ", nrow(bad), " library/libraries. ",
       "A low rate means the CellID convention is wrong, not that cells are ",
       "missing: cells the vendor annotated came from this same count step, so ",
       "they must be present. Check that colnames end in -1 and that the _<index> ",
       "suffix matches ARC_00's recovered ordering.")
}
ok(sprintf("every library matched at least 99%% of its loom cells (total %s)",
           format(sum(match_report$matched), big.mark = ",")))

gex_all <- do.call(cbind, gex_list)
rm(gex_list); invisible(gc(FALSE))
say("")
say(sprintf("  combined GEX: %s genes x %s cells",
            format(nrow(gex_all), big.mark = ","),
            format(ncol(gex_all), big.mark = ",")))

# Order the metadata to the matrix, so no row/column mismatch is possible.
setkey(md, CellID)
md_ord <- md[colnames(gex_all)]
stopifnot(identical(md_ord$CellID, colnames(gex_all)))

obj <- CreateSeuratObject(counts = gex_all, assay = "RNA", project = TISSUE)
rm(gex_all); invisible(gc(FALSE))

obj$sample      <- md_ord$Samples
obj$animal      <- base::sub("[NP]$", "", md_ord$Samples)
obj$timepoint   <- ifelse(grepl("P$", md_ord$Samples), "Pre",
                   ifelse(grepl("N$", md_ord$Samples), "Necropsy", "Single"))
obj$group       <- md_ord$Stimulation
obj$celltype    <- md_ord$celltype
obj$celltype_display <- md_ord$celltype_display
obj$library_index    <- md_ord$library_index
obj$tissue      <- TISSUE
# Arm A is THC/CBD, Arm B is VEH. Resolved by ARC_00 from the loom Stimulation
# attribute cross-tabulated against Samples.
obj$arm <- ifelse(grepl("THC", md_ord$Stimulation), "THC/CBD",
           ifelse(grepl("VEH", md_ord$Stimulation), "VEH", "Pre-infection"))

ok(paste0("Seurat object created with ", ncol(obj), " cells"))

subh("metadata as attached")
print(table(obj$group, obj$timepoint))
say("")
print(head(sort(table(obj$celltype_display), decreasing = TRUE), 10))

# The vendor's UMAP, so a rebuilt object can reproduce their existing figures
# rather than generating a new embedding.
umap <- as.matrix(md_ord[, .(umap_1, umap_2)])
rownames(umap) <- md_ord$CellID
colnames(umap) <- c("UMAP_1", "UMAP_2")
obj[["umap_vendor"]] <- CreateDimReducObject(
  embeddings = umap, key = "UMAPv_", assay = "RNA")
ok("vendor UMAP coordinates attached as reduction 'umap_vendor'")

# -----------------------------------------------------------------------------
# Cell 5: consensus peaks
# -----------------------------------------------------------------------------

hdr("Phase 3: consensus peak set")

peak_list <- list()
for (i in seq_along(libs)) {
  bed <- file.path(COUNT_ROOT, run_ids[i], "outs", "atac_peaks.bed")
  if (!file.exists(bed)) stop("not found: ", bed)
  d <- fread(bed, header = FALSE, select = 1:3,
             col.names = c("chrom", "start", "end"))
  # cellranger BED is 0-based half-open; GRanges is 1-based inclusive.
  gr <- GRanges(d$chrom, IRanges(d$start + 1L, d$end))
  peak_list[[i]] <- gr
  say(sprintf("  %-14s %8s peaks", run_ids[i], format(length(gr), big.mark = ",")))
}

pooled <- suppressWarnings(do.call(c, peak_list))
rm(peak_list); invisible(gc(FALSE))
say("")
say(sprintf("  pooled:          %9s", format(length(pooled), big.mark = ",")))

consensus <- reduce(pooled)
say(sprintf("  reduced union:   %9s", format(length(consensus), big.mark = ",")))

consensus <- consensus[as.character(seqnames(consensus)) %in% KEEP_CONTIGS]
say(sprintf("  standard contigs:%9s", format(length(consensus), big.mark = ",")))

w <- width(consensus)
consensus <- consensus[w >= PEAK_MIN_WIDTH & w <= PEAK_MAX_WIDTH]
say(sprintf("  width filtered:  %9s   (%d to %d bp)",
            format(length(consensus), big.mark = ","),
            PEAK_MIN_WIDTH, PEAK_MAX_WIDTH))

say("")
say(sprintf("  width: min %d, median %d, max %d",
            min(width(consensus)), as.integer(median(width(consensus))),
            max(width(consensus))))
say(sprintf("  contigs: %d", length(unique(as.character(seqnames(consensus))))))

# ARC_03 measured these against the vendor's universe. Restated as a check.
expected <- c(Brain = 204136, Colon = 345553)[TISSUE]
say("")
say(sprintf("  vendor's universe: %s   ratio %.3f",
            format(expected, big.mark = ","), length(consensus) / expected))
if (abs(length(consensus) / expected - 1) > 0.15) {
  wrn("the consensus set is more than 15% away from the vendor's universe, which ARC_03 did not see. Check the contig and width filters above.")
} else {
  ok("within 15% of the vendor's universe, matching what ARC_03 measured")
}

saveRDS(consensus, file.path(OUT_DIR, paste0("consensus_peaks_", TISSUE, ".rds")))
export(consensus, file.path(OUT_DIR, paste0("consensus_peaks_", TISSUE, ".bed")))
ok("consensus peaks written as .rds and .bed")

# -----------------------------------------------------------------------------
# Cell 6: fragment objects
#
# `cells` is a named vector: names are the object's cell names, values are the
# barcodes as they appear in the fragment file.
# -----------------------------------------------------------------------------

hdr("Phase 4: fragment objects")

frags <- vector("list", length(libs))
# The cells map is built from the object, not from the loom. Those two differ by
# the handful of loom cells our count run did not retain (one in Brain_40872),
# and a map naming a cell the object does not have produces a peak column with
# no expression column behind it.
obj_lib <- obj$library_index
names(obj_lib) <- colnames(obj)

for (i in seq_along(libs)) {
  fp <- file.path(COUNT_ROOT, run_ids[i], "outs", "atac_fragments.tsv.gz")
  if (!file.exists(fp)) stop("not found: ", fp)
  if (!file.exists(paste0(fp, ".tbi"))) wrn(paste0(run_ids[i], ": no .tbi index beside the fragments file"))

  obj_cells <- names(obj_lib)[obj_lib == i]
  file_bc <- base::sub("_\\d+$", "", obj_cells)   # strip the _<index> suffix
  cells_map <- setNames(file_bc, obj_cells)

  say(sprintf("  %-14s %6d cells  %s",
              run_ids[i], length(cells_map),
              format(structure(file.size(fp), class = "object_size"),
                     units = "auto")))

  frags[[i]] <- call_supported(
    CreateFragmentObject,
    required = list(path = fp, cells = cells_map),
    optional = list(validate.cells = VALIDATE_FRAGMENTS),
    label    = "CreateFragmentObject")
}
ok(paste0(length(frags), " fragment objects created"))

if (!"validate.cells" %in% names(formals(CreateFragmentObject))) {
  ok(paste0("Signac ", packageVersion("Signac"),
            " validates cells unconditionally, so VALIDATE_FRAGMENTS has no effect"))
} else if (VALIDATE_FRAGMENTS) {
  ok("cells validated against each fragment file")
} else {
  wrn("validate.cells was FALSE, so cell presence in the fragment files was not checked")
}

# -----------------------------------------------------------------------------
# Cell 7: FeatureMatrix, checkpointed
# -----------------------------------------------------------------------------

hdr("Phase 5: FeatureMatrix")

say(paste0("  quantifying ", format(length(consensus), big.mark = ","),
           " peaks across ", length(libs), " libraries"))
say("  The expensive step. Each library is checkpointed, so a timeout resumes.")

fm_list <- vector("list", length(libs))
for (i in seq_along(libs)) {
  ck <- file.path(CKPT_DIR, paste0(run_ids[i], "_featurematrix.rds"))
  if (USE_CHECKPOINTS && file.exists(ck)) {
    fm_list[[i]] <- readRDS(ck)
    ok(sprintf("%-14s from checkpoint (%d x %d)", run_ids[i],
               nrow(fm_list[[i]]), ncol(fm_list[[i]])))
    next
  }
  t0 <- Sys.time()
  # Same treatment as CreateFragmentObject: process_n and verbose are passed
  # only if this Signac version declares them.
  fm <- call_supported(
    FeatureMatrix,
    required = list(fragments = frags[[i]],
                    features  = consensus,
                    cells     = names(frags[[i]]@cells)),
    optional = list(process_n = FM_CHUNK, verbose = FALSE),
    label    = "FeatureMatrix")
  saveRDS(fm, ck)
  fm_list[[i]] <- fm
  say(sprintf("  %-14s %d x %d   %.1f min", run_ids[i], nrow(fm), ncol(fm),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(fm); invisible(gc(FALSE))
}

# Every library was quantified against the same consensus set, so the rows align.
rn <- rownames(fm_list[[1]])
for (i in seq_along(fm_list)) {
  if (!identical(rownames(fm_list[[i]]), rn)) {
    stop(run_ids[i], ": FeatureMatrix rownames differ from the first library. ",
         "A stale checkpoint from an earlier consensus set is the likely cause; ",
         "delete ", CKPT_DIR, " and rerun.")
  }
}

peak_mat <- do.call(cbind, fm_list)
rm(fm_list); invisible(gc(FALSE))
say("")
say(sprintf("  combined peak matrix: %s peaks x %s cells",
            format(nrow(peak_mat), big.mark = ","),
            format(ncol(peak_mat), big.mark = ",")))
say(sprintf("  nonzero entries: %s  (%.2f%% dense)",
            format(length(peak_mat@x), big.mark = ","),
            100 * length(peak_mat@x) / (as.numeric(nrow(peak_mat)) * ncol(peak_mat))))

common <- intersect(colnames(obj), colnames(peak_mat))
say("")
say(sprintf("  cells in both assays: %s of %s RNA cells",
            format(length(common), big.mark = ","),
            format(ncol(obj), big.mark = ",")))
if (length(common) < 0.99 * ncol(obj)) {
  wrn("more than 1% of RNA cells have no peak column. Subsetting to the intersection.")
}
obj <- subset(obj, cells = common)
peak_mat <- peak_mat[, colnames(obj), drop = FALSE]

# -----------------------------------------------------------------------------
# Cell 8: ChromatinAssay and gene annotation
# -----------------------------------------------------------------------------

hdr("Phase 6: ChromatinAssay and annotation")

say("  reading the GTF for the gene model track")
gtf <- rtracklayer::import(GTF)
gtf <- gtf[as.character(seqnames(gtf)) %in% KEEP_CONTIGS]
# Signac expects gene_name, gene_biotype, type and tx_id on the annotation.
if (!"gene_name" %in% names(mcols(gtf))) gtf$gene_name <- gtf$gene_id
gtf$gene_name <- ifelse(is.na(gtf$gene_name) | gtf$gene_name == "",
                        gtf$gene_id, gtf$gene_name)
if ("transcript_id" %in% names(mcols(gtf))) gtf$tx_id <- gtf$transcript_id
if (!"gene_biotype" %in% names(mcols(gtf))) gtf$gene_biotype <- "protein_coding"
say(sprintf("  annotation: %s records, %s genes",
            format(length(gtf), big.mark = ","),
            format(length(unique(gtf$gene_id)), big.mark = ",")))

chrom <- call_supported(
  CreateChromatinAssay,
  required = list(counts = peak_mat, fragments = frags, annotation = gtf),
  optional = list(genome = NULL, sep = c("-", "-")),
  label    = "CreateChromatinAssay")
obj[["ATAC"]] <- chrom
rm(peak_mat, chrom); invisible(gc(FALSE))
ok("ChromatinAssay added as assay 'ATAC'")

DefaultAssay(obj) <- "ATAC"
say("")
say(sprintf("  fragment objects attached: %d", length(Fragments(obj[["ATAC"]]))))
say(sprintf("  annotation records:        %s",
            format(length(Annotation(obj[["ATAC"]])), big.mark = ",")))

# Peak coordinates must use the vendor's chr1-56999-57879 style, or no join
# against their DA tables works.
say("")
say("  example peak names:")
say(paste0("    ", paste(head(rownames(obj[["ATAC"]]), 3), collapse = "  ")))
if (!grepl("^chr[0-9XY]+-[0-9]+-[0-9]+$", rownames(obj[["ATAC"]])[1])) {
  wrn("peak names are not in the vendor's chr-start-end form, so a join against their DA tables will fail")
} else {
  ok("peak names match the vendor's chr-start-end convention")
}

# -----------------------------------------------------------------------------
# Cell 9: per-cell ATAC QC
#
# ARC_03 flagged Brain_40875 at 1,783 median fragments per cell against 14,000
# to 20,000 elsewhere, and 40875 is a VEH animal in an arm of only four. These
# per-cell numbers are what a coverage track should be normalized against, so
# one thin library cannot dominate a group mean.
# -----------------------------------------------------------------------------

hdr("Phase 7: per-cell ATAC QC")

obj$atac_peak_counts <- colSums(GetAssayData(obj, assay = "ATAC", layer = "counts"))
obj$atac_peak_features <- colSums(GetAssayData(obj, assay = "ATAC", layer = "counts") > 0)

qc <- data.table(
  sample   = obj$sample,
  arm      = obj$arm,
  counts   = obj$atac_peak_counts,
  features = obj$atac_peak_features
)[, .(cells = .N,
      median_counts = as.integer(median(counts)),
      median_features = as.integer(median(features))),
  by = .(sample, arm)][order(sample)]

say("")
print(qc)

say("")
say("  Coverage tracks should be normalized PER ANIMAL rather than by pooling")
say("  fragments across a group. ARC_03 found one brain library an order of")
say("  magnitude below the others in fragments per cell, and with only four")
say("  animals per arm a single thin library is a quarter of that arm.")

low <- qc[median_counts < 0.4 * median(qc$median_counts)]
if (nrow(low)) {
  say("")
  for (i in seq_len(nrow(low))) {
    wrn(sprintf("%s (%s): median %d peak counts per cell, well below the tissue median of %d",
                low$sample[i], low$arm[i], low$median_counts[i],
                as.integer(median(qc$median_counts))))
  }
}

# -----------------------------------------------------------------------------
# Cell 10: save
# -----------------------------------------------------------------------------

hdr("Phase 8: save")

rds <- file.path(OUT_DIR, paste0("seurat_", TISSUE, "_rebuilt.rds"))
saveRDS(obj, rds)
ok(paste0("wrote ", rds, "  (",
          format(structure(file.size(rds), class = "object_size"), units = "auto"), ")"))

fwrite(match_report, file.path(OUT_DIR, paste0("barcode_match_", TISSUE, ".csv")))
fwrite(qc, file.path(OUT_DIR, paste0("atac_qc_", TISSUE, ".csv")))
ok("wrote barcode_match and atac_qc tables")

cm <- data.table(
  CellID = colnames(obj), sample = obj$sample, animal = obj$animal,
  timepoint = obj$timepoint, group = obj$group, arm = obj$arm,
  celltype = obj$celltype, celltype_display = obj$celltype_display,
  library_index = obj$library_index,
  atac_peak_counts = obj$atac_peak_counts,
  atac_peak_features = obj$atac_peak_features
)
fwrite(cm, file.path(OUT_DIR, paste0("cell_metadata_rebuilt_", TISSUE, ".csv")))
ok("wrote cell_metadata_rebuilt")

writeLines(capture.output({
  cat("ARC_04_build_object manifest\n")
  cat("tissue:            ", TISSUE, "\n")
  cat("built:             ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("host:              ", Sys.info()[["nodename"]], "\n")
  cat("elapsed (min):     ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "\n")
  cat("libraries:         ", paste(libs, collapse = ", "), "\n")
  cat("cells:             ", ncol(obj), "\n")
  cat("genes (RNA):       ", nrow(obj[["RNA"]]), "\n")
  cat("peaks (ATAC):      ", nrow(obj[["ATAC"]]), "\n")
  cat("fragment objects:  ", length(Fragments(obj[["ATAC"]])), "\n")
  cat("reference GTF:     ", GTF, "\n")
  cat("Seurat:            ", as.character(packageVersion("Seurat")), "\n")
  cat("Signac:            ", as.character(packageVersion("Signac")), "\n")
  cat("\nCell names follow the vendor's BARCODE-1_<index> convention, set\n")
  cat("explicitly rather than by merge(), so they join directly to the loom.\n")
  cat("Cells are an inner join against the loom annotation, so this object\n")
  cat("holds exactly the cells behind the collaborators' figures.\n")
}), file.path(OUT_DIR, paste0("MANIFEST_", TISSUE, ".txt")))
ok("wrote MANIFEST")

hdr("Done")
say(sprintf("  %s: %s cells, %s genes, %s peaks",
            TISSUE, format(ncol(obj), big.mark = ","),
            format(nrow(obj[["RNA"]]), big.mark = ","),
            format(nrow(obj[["ATAC"]]), big.mark = ",")))
say(sprintf("  elapsed: %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
say("")
say("  Ready for coverage tracks. CoveragePlot reads from the fragment objects,")
say("  so a gene window works directly:")
say("")
say("    obj <- readRDS(\"", rds, "\")")
say("    DefaultAssay(obj) <- \"ATAC\"")
say("    Idents(obj) <- \"celltype_display\"")
say("    CoveragePlot(obj, region = \"NRP1\", extend.upstream = 5000,")
say("                 extend.downstream = 5000, split.by = \"arm\")")
say("")
say("  Next: ARC_05_coverage_tracks.R for the gene panel the collaborators want.")

invisible(gc(FALSE))
