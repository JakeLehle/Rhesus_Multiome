#!/usr/bin/env Rscript
# =============================================================================
# ARC_05_coverage_tracks.R  <tissue> [gene ...]
#
# Genome-browser style ATAC coverage tracks for a named gene panel, restricted
# to one cell type and split by treatment arm. Reads the objects ARC_04 built.
#
#   Rscript ARC_05_coverage_tracks.R Brain
#   Rscript ARC_05_coverage_tracks.R Colon
#   Rscript ARC_05_coverage_tracks.R Brain NRP1          # one gene only
#
# On the cluster, the brain object is 1.2 GB and the colon object 1.7 GB, so:
#
#   sbatch --job-name=arc05_brain --partition=normal --cpus-per-task=4 \
#          --mem=96G --time=04:00:00 \
#          --output=/master/jlehle/WORKING/LOGS/arc05_brain_%j.out \
#          --wrap="source /master/jlehle/anaconda3/bin/activate multiome_R && \
#                  Rscript /master/jlehle/WORKING/SC/fastq/Rhesus_Multiome/ARC_05_coverage_tracks.R Brain"
#
# -----------------------------------------------------------------------------
# WHAT THIS PRODUCES PER GENE, AND WHY THREE PANELS RATHER THAN ONE
# -----------------------------------------------------------------------------
#
# 1. <GENE>_<tissue>_<celltype>_byArm         the deliverable. One track per
#    treatment arm across the gene body plus 20 kb upstream, with the gene
#    model and the consensus peak track underneath.
#
# 2. <GENE>_<tissue>_<celltype>_byAnimal      the honesty check. The same
#    window, one track per animal. CoveragePlot pools cells within a group and
#    normalizes by group depth, so an arm track is depth-weighted across its
#    animals. ARC_03 found Brain_40875 at 1,783 median fragments per cell
#    against 14,000 to 20,000 elsewhere, and 40875 is a VEH animal in an arm of
#    four. A real arm difference survives this panel. A difference carried by
#    one animal does not, and you want to know which before the figure leaves.
#
# 3. <GENE>_<tissue>_<celltype>_promoter      TSS plus or minus 2 kb, by arm.
#    Aim 1 asks whether a DEG is supported at the promoter or distally, and
#    that is unreadable at 25 kb width.
#
# Plus <GENE>_<tissue>_window_peaks.csv: every consensus peak in the window
# with its distance to the TSS and whether it sits in the promoter, so the
# figure can be described in text without re-measuring off the image.
#
# -----------------------------------------------------------------------------
# AND FOUR MORE PANELS WHERE A TIME COURSE EXISTS, WHICH IS COLON ONLY
# -----------------------------------------------------------------------------
#
# 4. _byTimepoint              pre-infection against post, pooled across arms
# 5. _byArmTimepoint           the four-way, control first: VEH before, VEH
#                              after, THC/CBD before, THC/CBD after
# 6. _byAnimalTimepoint        one pair of tracks per animal, pairs adjacent,
#                              control animals first
# 7. _promoter_byArmTimepoint  the four-way at TSS +/- 2 kb
#
# and, when PAIRED_ONLY leaves animals out, the same comparison over every
# animal that has cells at a timepoint:
#
# 8. _byArmTimepoint_allAnimals
# 9. _promoter_byArmTimepoint_allAnimals
#
# EVERY TRACK LABEL CARRIES ITS CELL COUNT, and the control arm is always the
# first track. A control arm that is thin is drawn and labelled n=6; a control
# arm that is genuinely absent produces a NO CONTROL TRACK warning in the log
# rather than a panel that looks complete.
#
# BRAIN HAS NO PRE-INFECTION LIBRARIES. All 8 are necropsy, which is forced by
# the tissue: you cannot biopsy basal ganglia from a living animal. The script
# says so and skips these four rather than drawing a single-track panel that
# reads like a result.
#
# TWO THINGS ABOUT THE COLON TIME COURSE THAT THE FIGURES CANNOT FIX.
#
# Only 4 of 8 colon animals have a pre-infection sample: 40874, 40875, 40877,
# 40878. By arm that is 3 THC/CBD against 1 VEH, so the VEH before-and-after
# rests on animal 40875 alone. A within-arm change there is one animal's
# biology and the figure cannot separate that from an arm effect. The script
# prints the pairing table and says this out loud every run.
#
# The timepoint panels use paired animals only by default (PAIRED_ONLY). With
# everything pooled, "pre" is 4 animals and "post" is 8, so part of any
# difference between the two tracks is just a difference in which animals are
# in them. Restricting to the paired set costs cells and buys the comparison.
#
# -----------------------------------------------------------------------------
# THE WINDOW, AND WHAT BOUNDS IT
# -----------------------------------------------------------------------------
#
# Default is 20 kb upstream and 5 kb downstream of the gene body. The 20 kb is
# not arbitrary: the vendor's own gene-to-peak map used a -20 kb to +2 kb
# window, so anything outside it has no counterpart in their DA tables. A
# distal element further out than that is visible here but unsupported by their
# statistics, which is worth saying out loud if one turns up.
#
# Jake Lehle, 2026-10
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Cell 1: the panel
# -----------------------------------------------------------------------------

ANALYSIS_ROOT <- "/master/jlehle/WORKING/SC/fastq/Rhesus_Multiome_ANALYSIS"
OBJ_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "04_objects")
OUT_ROOT      <- file.path(ANALYSIS_ROOT, "02_arc", "05_tracks")

# The vendor's tidy DA table from MULTIOME_01, used only to mark which peaks in
# the window they called differential. CHECK THIS PATH. If it is not there the
# script says so once and carries on without the annotation.
DA_PEAKS_CSV <- file.path(ANALYSIS_ROOT, "01_tidy", "peaks_da.csv")

# Cell type is matched as a regex against celltype_display, not compared for
# equality, because brain microglia are three clusters (Microglia-1 at 1,489
# cells, Microglia-2 at 124, Microglia-3 at 70). Pooling them gives roughly 840
# cells per arm; splitting leaves Microglia-3 at about 35 per arm, which is too
# thin to read as a coverage track. Per-cluster counts are printed so you can
# see what went in.
PANEL <- list(
  Brain = list(
    celltype_regex = "^Microglia",
    celltype_label = "Microglia",
    genes = c("NRP1", "SIGLEC9", "SIGLEC5", "AIM2")
  ),
  Colon = list(
    celltype_regex = "^Macrophage",
    celltype_label = "Macrophages",
    genes = c("HMOX1")
  )
)

# Symbol aliases tried in order. SIGLEG9 and SIGLEG5 as written in the request
# are kept as spellings to try, because a request is evidence about what the
# collaborator means, not about what Ensembl calls it.
GENE_ALIASES <- list(
  SIGLEC9 = c("SIGLEC9", "SIGLEG9", "SIGLEC-9"),
  SIGLEC5 = c("SIGLEC5", "SIGLEG5", "SIGLEC-5"),
  NRP1    = c("NRP1", "NRP-1", "NEUROPILIN1"),
  AIM2    = c("AIM2"),
  HMOX1   = c("HMOX1", "HO-1", "HMOX-1")
)
# When every alias misses, candidates are offered from this prefix.
GENE_FAMILY_HINT <- list(
  SIGLEC9 = "SIGLEC", SIGLEC5 = "SIGLEC",
  NRP1 = "NRP", AIM2 = "AIM", HMOX1 = "HMOX"
)

# Arms to draw, control first. Every panel orders tracks by this vector rather
# than alphabetically: "THC/CBD" sorts before "VEH" as a string, which would put
# treated above control and invert how the figure reads.
CONTROL_ARM <- "VEH"
ARMS_KEEP <- c("VEH", "THC/CBD")

# -----------------------------------------------------------------------------
# ARM COMES FROM THE ANIMAL, NOT FROM THE SAMPLE
#
# ARC_04 derived obj$arm from the loom's Stimulation attribute. At the
# pre-infection timepoint no treatment had been given, so those cells resolve to
# "Pre-infection" and carry no arm. That is correct as a description of the
# sample and useless for a paired before-and-after, where the question is what
# happened to THIS animal, which we know was later assigned to an arm.
#
# So arm is looked up per animal from the table below, and the lookup is
# cross-checked against the Stimulation of that animal's post-infection sample.
# A disagreement stops the run rather than mislabeling an arm.
#
# Source: ARC_00, from the loom Stimulation attribute cross-tabulated against
# Samples. Arm A is THC/CBD, Arm B is VEH.
# -----------------------------------------------------------------------------
ANIMAL_ARM <- c(
  "40872" = "THC/CBD", "40874" = "THC/CBD",
  "40877" = "THC/CBD", "40878" = "THC/CBD",
  "40873" = "VEH",     "40875" = "VEH",
  "40876" = "VEH",     "40879" = "VEH"
)

# ARC_04's timepoint values, and how they read on a figure. "Single" is brain,
# which has no pre-infection libraries at all.
TIMEPOINT_DISPLAY <- c(
  "Pre"      = "Pre-infection",
  "Necropsy" = "Post (necropsy)",
  "Single"   = "Single timepoint"
)
TIMEPOINT_ORDER <- c("Pre-infection", "Post (necropsy)")

# Restrict the timepoint figures to animals that have BOTH timepoints in this
# cell type. Otherwise "Pre" is four animals and "Post" is eight, and any
# difference between the two tracks is partly a difference in which animals are
# in them. Set FALSE to pool everything available and accept that confound.
PAIRED_ONLY <- TRUE

# -----------------------------------------------------------------------------
# THRESHOLDS, AND WHY THEY NO LONGER DELETE A TRACK GROUP
#
# An earlier version required 10 cells at both timepoints for an animal to count
# as paired. 40875 is the ONLY VEH animal with a pre-infection colon sample, so
# if its macrophages fell under that floor the whole VEH arm dropped out of the
# timepoint figure and the panel showed treated animals with nothing to compare
# them against. A threshold that can silently remove the control group is worse
# than no threshold.
#
# Now: an animal is paired if it has at least MIN_CELLS_PAIRED cells at both
# timepoints, which is deliberately 1. Thin groups are drawn, every track label
# carries its own n, and MIN_CELLS_FLAG decides only what gets called out in the
# log. The one remaining hard drop is per-animal tracks under MIN_CELLS_DROP,
# and any track it removes is named out loud.
# -----------------------------------------------------------------------------
MIN_CELLS_PAIRED <- 1
MIN_CELLS_FLAG   <- 50
MIN_CELLS_DROP   <- 5

# Track labels carry their cell count, so a 6-cell control track cannot be read
# as equivalent to a 600-cell one no matter who is looking at the figure.
SHOW_N_IN_LABELS <- TRUE

# With PAIRED_ONLY the arm-by-timepoint figure uses paired animals only, which
# is the clean within-animal comparison but can leave the control arm on one
# animal. This also emits an all-animals variant of that figure, which has more
# cells behind the control and a different confound: "pre" and "post" are then
# made of different animals. Two figures, two honest readings, named so they
# cannot be mixed up.
ALSO_ALL_ANIMALS <- TRUE

# Okabe-Ito, so the figures survive being printed in grayscale or read by
# someone with deuteranopia.
ARM_COLORS <- c("VEH" = "#0072B2", "THC/CBD" = "#D55E00",
                "Pre-infection" = "#999999")

# Timepoint uses value rather than hue: pre-infection is the lighter tone of the
# same family in both arms, so a reader tracks arm by color and timepoint by
# intensity instead of holding four unrelated hues in mind.
TP_COLORS <- c("Pre-infection" = "#999999", "Post (necropsy)" = "#333333")
ARM_TP_COLORS <- c(
  "VEH Pre-infection"       = "#9ECAE1",
  "VEH Post (necropsy)"     = "#0072B2",
  "THC/CBD Pre-infection"   = "#F0C9A8",
  "THC/CBD Post (necropsy)" = "#D55E00"
)

EXTEND_UP      <- 20000
EXTEND_DOWN    <- 5000
PROMOTER_FLANK <- 2000

FONT_SIZE <- 28          # project convention is 28 to 34
FIG_DPI   <- 300
FIG_W     <- 14
FIG_H_MIN <- 9

# Signac's default. Scales each group's coverage by its total reads in the
# assay, which is what makes the per-animal panel necessary rather than
# optional: it corrects depth between groups, not between animals pooled
# inside one group.
NORMALIZATION <- "RC"

# -----------------------------------------------------------------------------
# Cell 2: plumbing
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: ARC_05_coverage_tracks.R <Brain|Colon> [gene ...]")
TISSUE <- args[1]
if (!TISSUE %in% names(PANEL)) stop("tissue must be Brain or Colon")
GENES_REQ <- if (length(args) > 1) args[-1] else PANEL[[TISSUE]]$genes

say  <- function(...) cat(paste0(..., "\n"))
hdr  <- function(x) { say(""); say(strrep("=", 76)); say(x); say(strrep("=", 76)) }
subh <- function(x) { say(""); say(paste0("--- ", x, " ---")) }
ok   <- function(x) say(paste0("  [ok]    ", x))
wrn  <- function(x) say(paste0("  [warn]  ", x))

# Same guard as ARC_04, for the same reason: a helper named sub() shadows
# base::sub and the failure lands deep in a later phase.
for (fn in c("say", "hdr", "subh", "ok", "wrn")) {
  if (exists(fn, envir = baseenv(), inherits = FALSE)) {
    stop("helper '", fn, "' collides with base::", fn, ". Rename it.")
  }
}

# "VEH Pre-infection" plus its cell count. Applied to every track label in every
# panel so sample size travels with the picture.
lab_n <- function(base, n) {
  if (!SHOW_N_IN_LABELS) return(base)
  paste0(base, " (n=", format(n, big.mark = ","), ")")
}
strip_n <- function(x) base::sub(" \\(n=[0-9,]+\\)$", "", x)

# Colors are keyed by the undecorated label, so the palette has to be rebuilt
# against the decorated factor levels. If any level has no color, return NULL and
# let Signac's defaults stand rather than hand a manual scale an incomplete
# mapping, which fails at render time.
palette_for <- function(lvls, base_map) {
  base <- strip_n(lvls)
  out <- unname(base_map[base])
  if (anyNA(out)) return(NULL)
  names(out) <- lvls
  out
}

call_supported <- function(fn, required, optional = list(), label = "call") {
  f <- names(formals(fn))
  keep <- optional[names(optional) %in% f]
  dropped <- setdiff(names(optional), names(keep))
  if (length(dropped)) {
    wrn(paste0(label, ": this version does not take ",
               paste(dropped, collapse = ", "), "; its default applies"))
  }
  do.call(fn, c(required, keep))
}

CT_LABEL <- PANEL[[TISSUE]]$celltype_label
CT_SLUG  <- gsub("[^A-Za-z0-9]+", "", CT_LABEL)
OUT_DIR  <- file.path(OUT_ROOT, TISSUE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

t_start <- Sys.time()
hdr(paste0("ARC_05_coverage_tracks : ", TISSUE))
say(paste0("host:      ", Sys.info()[["nodename"]]))
say(paste0("started:   ", format(t_start, "%Y-%m-%d %H:%M:%S")))
say(paste0("cell type: ", CT_LABEL, "   (regex ", PANEL[[TISSUE]]$celltype_regex, ")"))
say(paste0("genes:     ", paste(GENES_REQ, collapse = ", ")))
say(paste0("window:    -", format(EXTEND_UP, big.mark = ","), " to +",
           format(EXTEND_DOWN, big.mark = ","), " bp around the gene body"))
say(paste0("out:       ", OUT_DIR))

# -----------------------------------------------------------------------------
# Cell 3: load, and confirm the fragment files are still reachable
#
# A ChromatinAssay stores a path and an MD5 per fragment file, not the data. If
# the count output moved since ARC_04 ran, every CoveragePlot would come back
# empty rather than erroring, so this is checked before anything is drawn.
# -----------------------------------------------------------------------------

hdr("Phase 1: object")

rds <- file.path(OBJ_ROOT, TISSUE, paste0("seurat_", TISSUE, "_rebuilt.rds"))
if (!file.exists(rds)) stop("not found: ", rds, "\nRun ARC_04_build_object.R first.")

say(paste0("  reading ", rds, " (",
           format(structure(file.size(rds), class = "object_size"), units = "auto"), ")"))
obj <- readRDS(rds)
DefaultAssay(obj) <- "ATAC"
ok(sprintf("%s cells, %s peaks, %s genes",
           format(ncol(obj), big.mark = ","),
           format(nrow(obj[["ATAC"]]), big.mark = ","),
           format(nrow(obj[["RNA"]]), big.mark = ",")))

frag_list <- Fragments(obj[["ATAC"]])
subh("fragment files")
missing_frag <- character(0)
for (f in frag_list) {
  p <- GetFragmentData(f, "path")
  n <- length(GetFragmentData(f, "cells"))
  if (file.exists(p)) {
    say(sprintf("  %6d cells  %s", n, p))
  } else {
    missing_frag <- c(missing_frag, p)
    wrn(sprintf("%6d cells  MISSING  %s", n, p))
  }
}
if (length(missing_frag)) {
  stop(length(missing_frag), " fragment file(s) are no longer at the path stored ",
       "in the object. Coverage plots would come back empty rather than failing, ",
       "so this stops here. Fix with Signac::UpdatePath() on each Fragment ",
       "object, or rerun ARC_04 if the count output moved for good.")
}
ok(paste0(length(frag_list), " fragment files reachable"))

ann <- Annotation(obj[["ATAC"]])
if (is.null(ann) || length(ann) == 0) stop("the ATAC assay carries no gene annotation")
ok(paste0(format(length(ann), big.mark = ","), " annotation records"))

# -----------------------------------------------------------------------------
# Cell 4: the cell group
# -----------------------------------------------------------------------------

hdr("Phase 2: cells")

ct <- as.character(obj$celltype_display)
hit <- grepl(PANEL[[TISSUE]]$celltype_regex, ct)
if (!any(hit)) {
  say("")
  say("  cell types present:")
  for (nm in sort(unique(ct))) say(paste0("    ", nm))
  stop("no cells matched '", PANEL[[TISSUE]]$celltype_regex, "'")
}

subh(paste0("clusters pooled into '", CT_LABEL, "'"))
for (nm in sort(unique(ct[hit]))) {
  say(sprintf("  %-32s %6d cells", nm, sum(ct == nm)))
}
ok(sprintf("%d cells matched", sum(hit)))

sub_obj <- subset(obj, cells = colnames(obj)[hit])
rm(obj); invisible(gc(FALSE))

# --- arm, from the animal ----------------------------------------------------

animal <- as.character(sub_obj$animal)
unknown <- setdiff(unique(animal), names(ANIMAL_ARM))
if (length(unknown)) {
  stop("animal(s) not in ANIMAL_ARM: ", paste(unknown, collapse = ", "),
       ". Add them before running, or the arm labels are guesses.")
}
arm_animal <- unname(ANIMAL_ARM[animal])

# Cross-check against the loom's own Stimulation, using post-infection samples
# only, since a pre-infection sample has no treatment to disagree about.
subh("arm lookup cross-check")
tp_raw <- as.character(sub_obj$timepoint)
chk <- data.table(animal = animal, stim = as.character(sub_obj$group),
                  tp = tp_raw, arm = arm_animal)[tp != "Pre"]
bad_arm <- character(0)
for (a in sort(unique(chk$animal))) {
  stims <- unique(chk[animal == a, stim])
  expect <- ANIMAL_ARM[[a]]
  implied <- if (any(grepl("THC", stims))) "THC/CBD" else
             if (any(grepl("VEH", stims))) "VEH" else NA_character_
  flag <- if (is.na(implied)) "no THC/VEH in Stimulation" else
          if (implied != expect) "MISMATCH" else "ok"
  if (flag == "MISMATCH") bad_arm <- c(bad_arm, a)
  say(sprintf("  %-8s table says %-8s  loom says %-24s  %s",
              a, expect, paste(stims, collapse = "/"), flag))
}
if (length(bad_arm)) {
  stop("ANIMAL_ARM disagrees with the loom Stimulation for: ",
       paste(bad_arm, collapse = ", "),
       ". Resolve this before any figure goes out; an arm label is the whole ",
       "claim of these panels.")
}
ok("ANIMAL_ARM agrees with the loom for every post-infection sample")

arms_present <- intersect(ARMS_KEEP, unique(arm_animal))
if (length(arms_present) < 2) {
  wrn(paste0("only ", length(arms_present), " of the requested arms are present"))
}
keep_cells <- colnames(sub_obj)[arm_animal %in% arms_present]
sub_obj <- subset(sub_obj, cells = keep_cells)

sub_obj$arm <- factor(unname(ANIMAL_ARM[as.character(sub_obj$animal)]),
                      levels = arms_present)

# --- timepoint ---------------------------------------------------------------

tp_raw <- as.character(sub_obj$timepoint)
unknown_tp <- setdiff(unique(tp_raw), names(TIMEPOINT_DISPLAY))
if (length(unknown_tp)) {
  stop("unrecognized timepoint value(s): ", paste(unknown_tp, collapse = ", "))
}
tp_disp <- unname(TIMEPOINT_DISPLAY[tp_raw])
tp_levels <- intersect(TIMEPOINT_ORDER, unique(tp_disp))
sub_obj$tp <- factor(tp_disp, levels = c(tp_levels, setdiff(unique(tp_disp), tp_levels)))

subh("cells per timepoint")
for (nm in levels(sub_obj$tp)) {
  if (sum(tp_disp == nm) > 0) say(sprintf("  %-20s %6d cells", nm, sum(tp_disp == nm)))
}

HAS_TIMECOURSE <- length(tp_levels) >= 2

# Cells for the three panels that were already approved. Those compare arms at
# one timepoint, so pre-infection cells stay out of them. For brain this is
# every cell; for colon it is the necropsy samples, which is exactly the set the
# earlier version happened to use when pre-infection cells had no arm.
post_cells <- colnames(sub_obj)[as.character(sub_obj$tp) != "Pre-infection"]
POST_OBJ <- subset(sub_obj, cells = post_cells)
POST_OBJ$arm <- droplevels(POST_OBJ$arm)
arms_post <- levels(POST_OBJ$arm)

# Arm tracks carry their n. The arm factor was built with levels in ARMS_KEEP
# order, so control comes first.
arm_counts <- table(POST_OBJ$arm)
POST_OBJ$arm_lab <- factor(
  lab_n(as.character(POST_OBJ$arm), arm_counts[as.character(POST_OBJ$arm)]),
  levels = lab_n(levels(POST_OBJ$arm), arm_counts[levels(POST_OBJ$arm)]))

# Ordered by the FACTOR, not by as.character: "THC/CBD" sorts before "VEH" as a
# string, which would stack treated animals above controls.
at_base <- paste0(POST_OBJ$sample, " ", as.character(POST_OBJ$arm))
at_cnt <- table(at_base)
POST_OBJ$animal_track <- factor(
  lab_n(at_base, at_cnt[at_base]),
  levels = lab_n(unique(at_base[order(POST_OBJ$arm, POST_OBJ$sample)]),
                 at_cnt[unique(at_base[order(POST_OBJ$arm, POST_OBJ$sample)])]))

subh("cells per animal track (arm panels)")
at <- table(POST_OBJ$animal_track)
for (nm in names(at)) say(sprintf("  %-34s %6d cells", nm, at[[nm]]))
thin <- names(at)[at < MIN_CELLS_FLAG]
if (length(thin)) {
  wrn(paste0(length(thin), " track(s) under ", MIN_CELLS_FLAG, " cells: ",
             paste(strip_n(thin), collapse = ", "), ". Indicative only."))
}
say("")
say(paste0("  arm track order: ", paste(levels(POST_OBJ$arm_lab), collapse = "  then  ")))
if (!CONTROL_ARM %in% levels(POST_OBJ$arm)) {
  wrn(paste0("the control arm '", CONTROL_ARM, "' is absent from this panel, so ",
             "there is nothing to compare the treated arm against"))
} else if (levels(POST_OBJ$arm)[1] != CONTROL_ARM) {
  wrn(paste0("control arm '", CONTROL_ARM, "' is not first; check ARMS_KEEP"))
} else {
  ok(paste0("control arm '", CONTROL_ARM, "' is first"))
}

# --- the paired timepoint set ------------------------------------------------

TP_OBJ <- NULL
TP_ALL_OBJ <- NULL
paired_animals <- character(0)
if (!HAS_TIMECOURSE) {
  subh("timepoint panels")
  wrn(paste0(TISSUE, " has a single timepoint (", paste(tp_levels, collapse = ", "),
             "), so there is no before-and-after to draw."))
  if (TISSUE == "Brain") {
    say("          All 8 brain libraries are necropsy-only: the vendor collected")
    say("          no pre-infection brain, which is anatomically unavoidable in a")
    say("          living animal. The timepoint figures are colon only, and no")
    say("          amount of regrouping creates a brain pre-infection track.")
    say("")
    say("          The control comparison still exists in brain: _byArm puts")
    say("          VEH (4 animals) first and THC/CBD (4 animals) second, which")
    say("          is a balanced control-versus-treated contrast. What brain")
    say("          cannot have is the time axis.")
  }
} else {
  tab_at <- table(as.character(sub_obj$animal), as.character(sub_obj$tp))
  present <- tab_at[, intersect(tp_levels, colnames(tab_at)), drop = FALSE] >= MIN_CELLS_PAIRED
  paired_animals <- rownames(present)[rowSums(present) == length(tp_levels)]

  subh("timepoint coverage per animal")
  say(sprintf("  %-8s %-6s %12s %16s %s", "animal", "arm", "pre", "post", "paired"))
  for (a in sort(rownames(tab_at))) {
    pre <- if ("Pre-infection" %in% colnames(tab_at)) tab_at[a, "Pre-infection"] else 0
    pos <- if ("Post (necropsy)" %in% colnames(tab_at)) tab_at[a, "Post (necropsy)"] else 0
    say(sprintf("  %-8s %-8s %10d %16d %s", a, ANIMAL_ARM[[a]], pre, pos,
                if (a %in% paired_animals) "yes" else "no"))
  }

  if (!length(paired_animals)) {
    wrn(paste0("no animal has at least ", MIN_CELLS_DROP, " cells at both ",
               "timepoints in this cell type, so the timepoint panels are skipped"))
  } else {
    by_arm <- table(unname(ANIMAL_ARM[paired_animals]))
    say("")
        ok(paste0(length(paired_animals), " paired animal(s): ",
              paste(paired_animals, collapse = ", ")))
    say(paste0("          by arm: ",
               paste(sprintf("%s n=%d", names(by_arm), as.integer(by_arm)),
                     collapse = ", ")))

    # The imbalance is the headline caveat on these panels, so it is stated
    # here rather than left for someone to work out from the table above.
    if (any(by_arm < 2)) {
      thin_arm <- names(by_arm)[by_arm < 2]
      wrn(paste0("the ", paste(thin_arm, collapse = " and "),
                 " before-and-after rests on a single animal."))
      say("          A within-arm change there is one animal's biology, not an")
      say("          arm effect, and the figure cannot distinguish the two.")
      say("          Worth raising with the collaborators before this is shown")
      say("          as a treatment result.")
    }
    if (length(by_arm) < length(arms_present)) {
      wrn(paste0("no paired animals at all in: ",
                 paste(setdiff(arms_present, names(by_arm)), collapse = ", ")))
    }

    # Builds the arm-by-timepoint and animal-by-timepoint factors on a subset.
    # Levels are enumerated from ARMS_KEEP and the timepoint order rather than
    # from what happens to be in the data, and droplevels is applied ONLY to
    # levels with zero cells. That is the difference between a control arm with
    # six cells, which gets drawn and labelled n=6, and a control arm that
    # quietly disappears from the figure.
    build_tp <- function(cells, animals_for_order, tag) {
      o <- subset(sub_obj, cells = cells)
      o$tp  <- droplevels(factor(as.character(o$tp), levels = levels(sub_obj$tp)))
      o$arm <- factor(as.character(o$arm), levels = arms_present)

      c_tp <- table(as.character(o$tp))
      o$tp_lab <- factor(lab_n(as.character(o$tp), c_tp[as.character(o$tp)]),
                         levels = lab_n(levels(o$tp), c_tp[levels(o$tp)]))

      base_at <- paste(as.character(o$arm), as.character(o$tp))
      lv_at <- unlist(lapply(arms_present, function(a) paste(a, levels(o$tp))))
      lv_at <- lv_at[lv_at %in% base_at]          # zero-cell levels only
      c_at <- table(base_at)
      o$arm_tp <- factor(lab_n(base_at, c_at[base_at]),
                         levels = lab_n(lv_at, c_at[lv_at]))

      # Animal first, then timepoint, so each animal's pair sits adjacent and
      # the eye reads before-to-after down the panel. Animals ordered by arm
      # using the factor, so control animals come first.
      ord <- animals_for_order[order(
        factor(unname(ANIMAL_ARM[animals_for_order]), levels = arms_present),
        animals_for_order)]
      base_an <- paste0(as.character(o$animal), " ", as.character(o$tp))
      lv_an <- unlist(lapply(ord, function(a) paste0(a, " ", levels(o$tp))))
      lv_an <- lv_an[lv_an %in% base_an]
      c_an <- table(base_an)

      drop_an <- names(c_an)[c_an < MIN_CELLS_DROP]
      if (length(drop_an)) {
        wrn(paste0(tag, ": dropping ", length(drop_an),
                   " animal track(s) under ", MIN_CELLS_DROP, " cells: ",
                   paste(sprintf("%s n=%d", drop_an, as.integer(c_an[drop_an])),
                         collapse = ", ")))
        lv_an <- setdiff(lv_an, drop_an)
      }
      o$animal_tp <- factor(lab_n(base_an, c_an[base_an]),
                            levels = lab_n(lv_an, c_an[lv_an]))
      o
    }

    report_tp <- function(o, tag) {
      subh(paste0("cells per track: ", tag))
      for (nm in levels(o$arm_tp)) {
        n <- sum(as.character(o$arm_tp) == nm)
        say(sprintf("  %-40s %6d cells%s", nm, n,
                    if (n < MIN_CELLS_FLAG) "   (thin)" else ""))
      }
      present_arms <- unique(strip_n(sub(" (Pre-infection|Post \\(necropsy\\))$", "",
                                         levels(o$arm_tp))))
      if (!CONTROL_ARM %in% present_arms) {
        wrn(paste0("NO CONTROL TRACK: '", CONTROL_ARM, "' has no cells at either ",
                   "timepoint in this cell type, so this panel shows the treated ",
                   "arm with nothing to compare it against. Do not read a ",
                   "before-and-after change off it."))
      } else {
        ok(paste0("control arm '", CONTROL_ARM, "' present and ordered first"))
      }
      say("")
      for (nm in levels(o$animal_tp)) {
        n <- sum(as.character(o$animal_tp) == nm)
        say(sprintf("  %-40s %6d cells%s", nm, n,
                    if (n < MIN_CELLS_FLAG) "   (thin)" else ""))
      }
    }

    tp_cells <- if (PAIRED_ONLY) {
      colnames(sub_obj)[as.character(sub_obj$animal) %in% paired_animals]
    } else {
      colnames(sub_obj)
    }
    TP_OBJ <- build_tp(tp_cells, paired_animals, "paired animals")
    report_tp(TP_OBJ, "paired animals")
    ok(sprintf("%s cells across %d timepoint tracks and %d animal tracks",
               format(ncol(TP_OBJ), big.mark = ","),
               length(levels(TP_OBJ$arm_tp)), length(levels(TP_OBJ$animal_tp))))

    # The all-animals variant. More cells behind the control, at the cost of
    # "pre" and "post" being made of different animals. Only built when it would
    # actually differ from the paired version.
    if (ALSO_ALL_ANIMALS && PAIRED_ONLY &&
        length(paired_animals) < length(unique(as.character(sub_obj$animal)))) {
      TP_ALL_OBJ <- build_tp(colnames(sub_obj),
                             unique(as.character(sub_obj$animal)), "all animals")
      report_tp(TP_ALL_OBJ, "all animals")
      wrn("in the all-animals panel, 'pre' and 'post' are different animal sets,")
      say("          so part of any difference between the two is which animals")
      say("          are in each track. Use it for the control comparison, and")
      say("          the paired panel for the within-animal claim.")
    }
  }
}

ok(sprintf("%s cells in %d arms for the arm panels",
           format(ncol(POST_OBJ), big.mark = ","), length(arms_post)))

# -----------------------------------------------------------------------------
# Cell 5: resolve gene symbols against the annotation we actually have
#
# Ensembl Mmul_8.0.1 carries plenty of non-standard symbols, so a miss here is
# expected rather than exceptional. Candidates are printed from the family
# prefix so the next run can name the gene correctly instead of guessing.
# -----------------------------------------------------------------------------

hdr("Phase 3: gene symbols")

ann_names <- unique(as.character(ann$gene_name))
ann_names <- ann_names[!is.na(ann_names) & ann_names != ""]

resolve_gene <- function(g) {
  tries <- unique(c(g, GENE_ALIASES[[g]]))
  for (t in tries) {
    if (t %in% ann_names) return(list(found = TRUE, symbol = t, tried = tries))
  }
  # case-insensitive second pass, for symbols that differ only in case
  for (t in tries) {
    i <- which(toupper(ann_names) == toupper(t))
    if (length(i)) return(list(found = TRUE, symbol = ann_names[i[1]], tried = tries))
  }
  list(found = FALSE, symbol = NA_character_, tried = tries)
}

resolved <- list()
for (g in GENES_REQ) {
  r <- resolve_gene(g)
  if (r$found) {
    if (r$symbol != g) {
      ok(sprintf("%-10s found in the annotation as '%s'", g, r$symbol))
    } else {
      ok(sprintf("%-10s found", g))
    }
    resolved[[g]] <- r$symbol
  } else {
    wrn(sprintf("%-10s not in the annotation. Tried: %s", g,
                paste(r$tried, collapse = ", ")))
    hint <- GENE_FAMILY_HINT[[g]]
    if (!is.null(hint)) {
      cand <- sort(grep(paste0("^", hint), ann_names, value = TRUE,
                        ignore.case = TRUE))
      if (length(cand)) {
        say(paste0("            candidates starting '", hint, "': ",
                   paste(head(cand, 25), collapse = ", "),
                   if (length(cand) > 25) " ..." else ""))
      } else {
        say(paste0("            nothing in the annotation starts '", hint, "'"))
      }
    }
    # Is it in the RNA assay but not the ATAC annotation? Different problem.
    if (g %in% rownames(sub_obj[["RNA"]])) {
      say("            note: this symbol IS a feature in the RNA assay, so the")
      say("            gap is in the GTF annotation attached to the ATAC assay,")
      say("            not in the reference as a whole.")
    }
  }
}

if (!length(resolved)) {
  stop("none of the requested genes resolved against the annotation. The ",
       "candidate lists above are the place to start.")
}
ok(paste0(length(resolved), " of ", length(GENES_REQ), " genes resolved"))

# -----------------------------------------------------------------------------
# Cell 6: the vendor's DA peaks, if the tidy table is there
# -----------------------------------------------------------------------------

da_gr <- NULL
if (file.exists(DA_PEAKS_CSV)) {
  da_ok <- tryCatch({
    da <- fread(DA_PEAKS_CSV)
    pk_col <- intersect(c("atac_peak", "peak", "region"), names(da))[1]
    if (is.na(pk_col)) stop("no peak column")
    p <- as.character(da[[pk_col]])
    m <- regmatches(p, regexec("^(.*)-(\\d+)-(\\d+)$", p))
    good <- vapply(m, length, integer(1)) == 4
    gr <- GRanges(vapply(m[good], `[`, character(1), 2),
                  IRanges(as.numeric(vapply(m[good], `[`, character(1), 3)),
                          as.numeric(vapply(m[good], `[`, character(1), 4))))
    ct_col <- intersect(c("atac_celltype", "celltype", "cell_type"), names(da))[1]
    if (!is.na(ct_col)) gr$celltype <- as.character(da[[ct_col]])[good]
    gr
  }, error = function(e) { wrn(paste0("could not use ", DA_PEAKS_CSV, ": ", conditionMessage(e))); NULL })
  if (!is.null(da_ok)) {
    da_gr <- da_ok
    ok(paste0(format(length(da_gr), big.mark = ","),
              " vendor DA peaks loaded for window annotation"))
  }
} else {
  wrn(paste0("no DA peak table at ", DA_PEAKS_CSV,
             "; window peak tables will not be marked for differential calls"))
}

# -----------------------------------------------------------------------------
# Cell 7: draw
# -----------------------------------------------------------------------------

hdr("Phase 4: coverage tracks")

Idents(POST_OBJ) <- "arm"
if (!is.null(TP_OBJ)) Idents(TP_OBJ) <- "arm_tp"

gene_region <- function(symbol) {
  g <- ann[!is.na(ann$gene_name) & ann$gene_name == symbol]
  if (!length(g)) return(NULL)
  # One gene can appear as many transcript and exon records; take the span.
  gr <- GRanges(as.character(seqnames(g))[1],
                IRanges(min(start(g)), max(end(g))),
                strand = as.character(strand(g))[1])
  gr
}

tss_of <- function(gr) {
  if (as.character(strand(gr)) == "-") end(gr) else start(gr)
}

save_fig <- function(p, stem, n_tracks) {
  h <- max(FIG_H_MIN, 3 + 1.6 * n_tracks)
  pdf_path <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(stem, ".png"))
  ggsave(pdf_path, p, width = FIG_W, height = h, units = "in",
         device = grDevices::cairo_pdf, limitsize = FALSE)
  ggsave(png_path, p, width = FIG_W, height = h, units = "in",
         dpi = FIG_DPI, limitsize = FALSE)
  say(sprintf("    wrote %s.pdf and .png  (%.0f x %.0f in, %d dpi)",
              stem, FIG_W, h, FIG_DPI))
}

# Gene names are italic, everything else roman. The gene goes in the title on
# its own so the whole title can be italicized; the cell type and grouping go
# in the subtitle, which stays roman.
decorate <- function(p, gene, subtitle) {
  p +
    plot_annotation(
      title = gene,
      subtitle = subtitle,
      theme = theme(
        plot.title    = element_text(size = FONT_SIZE + 4, face = "bold.italic"),
        plot.subtitle = element_text(size = FONT_SIZE)
      )
    ) &
    theme(text = element_text(size = FONT_SIZE))
}

# Project convention is explicit hex codes, but whether CoveragePlot accepts a
# cols argument depends on the Signac version, and a manual fill scale that does
# not match the data fails at render time rather than when it is added, which
# would cost the figure instead of just the colors. So: pass cols through the
# formals filter, and if this version ignored it, try the fill scale and PROVE
# it renders on a throwaway copy before keeping it. A recolor that cannot render
# is dropped and Signac's defaults stand.
apply_arm_colors <- function(p, cols) {
  if (is.null(cols) || !length(cols)) return(p)
  cand <- tryCatch(p & scale_fill_manual(values = cols, drop = TRUE),
                   error = function(e) NULL)
  if (is.null(cand)) return(p)
  tmp <- tempfile(fileext = ".png")
  rendered <- tryCatch({
    suppressWarnings(ggsave(tmp, cand, width = 4, height = 3, dpi = 36))
    TRUE
  }, error = function(e) FALSE)
  unlink(tmp)
  if (rendered) cand else {
    wrn("the manual arm palette would not render; keeping Signac's default colors")
    p
  }
}

draw <- function(object, region, group_by, cols, extend_up, extend_down) {
  p <- call_supported(
    CoveragePlot,
    required = list(object = object, region = region),
    optional = list(
      group.by          = group_by,
      cols              = cols,
      extend.upstream   = extend_up,
      extend.downstream = extend_down,
      annotation        = TRUE,
      peaks             = TRUE,
      tile              = FALSE,
      links             = FALSE,       # no Links in this object; TRUE warns
      normalization     = NORMALIZATION
    ),
    label = "CoveragePlot")
  if (!"cols" %in% names(formals(CoveragePlot))) p <- apply_arm_colors(p, cols)
  p
}

window_peaks_table <- function(symbol, region, extend_up, extend_down) {
  w <- region
  if (as.character(strand(w)) == "-") {
    start(w) <- max(1, start(w) - extend_down); end(w) <- end(w) + extend_up
  } else {
    start(w) <- max(1, start(w) - extend_up);   end(w) <- end(w) + extend_down
  }
  pk <- StringToGRanges(rownames(POST_OBJ[["ATAC"]]), sep = c("-", "-"))
  idx <- which(as.logical(countOverlaps(pk, w) > 0))
  if (!length(idx)) return(NULL)
  tss <- tss_of(region)
  mid <- (start(pk[idx]) + end(pk[idx])) %/% 2
  dist <- mid - tss
  if (as.character(strand(region)) == "-") dist <- -dist
  dt <- data.table(
    gene = symbol,
    peak = rownames(POST_OBJ[["ATAC"]])[idx],
    chrom = as.character(seqnames(pk[idx])),
    start = start(pk[idx]), end = end(pk[idx]),
    width = width(pk[idx]),
    dist_to_tss = dist,
    region_class = ifelse(abs(dist) <= PROMOTER_FLANK, "promoter",
                   ifelse(dist < 0, "upstream_distal", "downstream_distal"))
  )
  if (!is.null(da_gr)) {
    dt$vendor_da_peak <- countOverlaps(pk[idx], da_gr) > 0
  }
  dt[order(dist_to_tss)]
}

made <- data.table()
for (g in names(resolved)) {
  symbol <- resolved[[g]]
  region <- gene_region(symbol)
  subh(paste0(g, "  (", symbol, ")"))

  if (is.null(region)) { wrn("could not build a region; skipped"); next }
  say(sprintf("  %s:%s-%s  (%s strand, %s bp)",
              as.character(seqnames(region)),
              format(start(region), big.mark = ","),
              format(end(region), big.mark = ","),
              as.character(strand(region)),
              format(width(region), big.mark = ",")))

  stem_base <- paste0(g, "_", TISSUE, "_", CT_SLUG)

  prom <- region
  t0 <- tss_of(region)
  start(prom) <- max(1, t0 - PROMOTER_FLANK); end(prom) <- t0 + PROMOTER_FLANK
  strand(prom) <- strand(region)

  tp_note <- if (PAIRED_ONLY) "paired animals only" else "all available animals"

  res <- tryCatch({
    # --- the arm panels, post-infection cells only ---------------------------
    # Control track first in all three, because the figure is read as treated
    # against control and not the other way round.
    # 1. by arm, wide
    p1 <- draw(POST_OBJ, region, "arm_lab",
               palette_for(levels(POST_OBJ$arm_lab), ARM_COLORS),
               EXTEND_UP, EXTEND_DOWN)
    p1 <- decorate(p1, g, sprintf("%s %s, control then treated, %s normalization",
                                  TISSUE, CT_LABEL, NORMALIZATION))
    save_fig(p1, paste0(stem_base, "_byArm"), length(levels(POST_OBJ$arm_lab)))

    # 2. by animal, same window
    p2 <- draw(POST_OBJ, region, "animal_track", NULL, EXTEND_UP, EXTEND_DOWN)
    p2 <- decorate(p2, g, sprintf("%s %s, one track per animal, control animals first",
                                  TISSUE, CT_LABEL))
    save_fig(p2, paste0(stem_base, "_byAnimal"), length(at))

    # 3. promoter zoom, by arm
    p3 <- draw(POST_OBJ, prom, "arm_lab",
               palette_for(levels(POST_OBJ$arm_lab), ARM_COLORS), 0, 0)
    p3 <- decorate(p3, g, sprintf("%s %s, TSS +/- %s bp, control then treated",
                                  TISSUE, CT_LABEL,
                                  format(PROMOTER_FLANK, big.mark = ",")))
    save_fig(p3, paste0(stem_base, "_promoter"), length(levels(POST_OBJ$arm_lab)))

    # --- the timepoint panels, colon only ------------------------------------
    if (!is.null(TP_OBJ)) {
      # 4. pre vs post, pooled across arms. The plainest version of the
      # question the collaborators asked, and the one with the most cells
      # behind each track.
      p4 <- draw(TP_OBJ, region, "tp_lab",
                 palette_for(levels(TP_OBJ$tp_lab), TP_COLORS),
                 EXTEND_UP, EXTEND_DOWN)
      p4 <- decorate(p4, g, sprintf("%s %s, pre vs post infection, %s",
                                    TISSUE, CT_LABEL, tp_note))
      save_fig(p4, paste0(stem_base, "_byTimepoint"), length(levels(TP_OBJ$tp)))

      # 5. arm by timepoint. Four tracks in control-first order: VEH before,
      # VEH after, THC/CBD before, THC/CBD after. This is the panel the
      # collaborators asked for, and the one where the 3:1 split in which
      # animals have a pre-infection sample bites hardest.
      p5 <- draw(TP_OBJ, region, "arm_tp",
                 palette_for(levels(TP_OBJ$arm_tp), ARM_TP_COLORS),
                 EXTEND_UP, EXTEND_DOWN)
      p5 <- decorate(p5, g, sprintf("%s %s, control then treated, before and after, %s",
                                    TISSUE, CT_LABEL, tp_note))
      save_fig(p5, paste0(stem_base, "_byArmTimepoint"),
               length(levels(TP_OBJ$arm_tp)))

      # 6. per animal, pairs adjacent. The only panel in which a within-animal
      # change is visible as such, rather than inferred from two pooled tracks
      # made of different cells.
      p6 <- draw(TP_OBJ, region, "animal_tp", NULL, EXTEND_UP, EXTEND_DOWN)
      p6 <- decorate(p6, g, sprintf("%s %s, before and after within each animal",
                                    TISSUE, CT_LABEL))
      save_fig(p6, paste0(stem_base, "_byAnimalTimepoint"),
               length(levels(TP_OBJ$animal_tp)))

      # 7. promoter zoom, arm by timepoint
      p7 <- draw(TP_OBJ, prom, "arm_tp",
                 palette_for(levels(TP_OBJ$arm_tp), ARM_TP_COLORS), 0, 0)
      p7 <- decorate(p7, g, sprintf("%s %s, TSS +/- %s bp, control then treated, before and after",
                                    TISSUE, CT_LABEL,
                                    format(PROMOTER_FLANK, big.mark = ",")))
      save_fig(p7, paste0(stem_base, "_promoter_byArmTimepoint"),
               length(levels(TP_OBJ$arm_tp)))
    }

    # 8 and 9. The all-animals variant of the arm-by-timepoint comparison, with
    # more cells behind the control arm and the different-animal-sets confound
    # in its place. Named so it cannot be confused with the paired version.
    if (!is.null(TP_ALL_OBJ)) {
      p8 <- draw(TP_ALL_OBJ, region, "arm_tp",
                 palette_for(levels(TP_ALL_OBJ$arm_tp), ARM_TP_COLORS),
                 EXTEND_UP, EXTEND_DOWN)
      p8 <- decorate(p8, g, sprintf("%s %s, control then treated, before and after, all animals",
                                    TISSUE, CT_LABEL))
      save_fig(p8, paste0(stem_base, "_byArmTimepoint_allAnimals"),
               length(levels(TP_ALL_OBJ$arm_tp)))

      p9 <- draw(TP_ALL_OBJ, prom, "arm_tp",
                 palette_for(levels(TP_ALL_OBJ$arm_tp), ARM_TP_COLORS), 0, 0)
      p9 <- decorate(p9, g, sprintf("%s %s, TSS +/- %s bp, before and after, all animals",
                                    TISSUE, CT_LABEL,
                                    format(PROMOTER_FLANK, big.mark = ",")))
      save_fig(p9, paste0(stem_base, "_promoter_byArmTimepoint_allAnimals"),
               length(levels(TP_ALL_OBJ$arm_tp)))
    }

    # peak table for the window
    wt <- window_peaks_table(symbol, region, EXTEND_UP, EXTEND_DOWN)
    if (!is.null(wt)) {
      fwrite(wt, file.path(OUT_DIR, paste0(g, "_", TISSUE, "_window_peaks.csv")))
      say(sprintf("    %d consensus peaks in the window: %d promoter, %d upstream distal, %d downstream distal",
                  nrow(wt), sum(wt$region_class == "promoter"),
                  sum(wt$region_class == "upstream_distal"),
                  sum(wt$region_class == "downstream_distal")))
      if (!is.null(wt$vendor_da_peak)) {
        say(sprintf("    %d of them overlap a vendor DA peak", sum(wt$vendor_da_peak)))
      }
    } else {
      wrn("no consensus peaks in this window, so the tracks will be flat")
    }
    TRUE
  }, error = function(e) { wrn(paste0("failed: ", conditionMessage(e))); FALSE })

  made <- rbind(made, data.table(gene = g, symbol = symbol,
                                 region = paste0(as.character(seqnames(region)), ":",
                                                 start(region), "-", end(region)),
                                 ok = res))
}

# -----------------------------------------------------------------------------
# Cell 8: summary
# -----------------------------------------------------------------------------

hdr("Done")
print(made)
fwrite(made, file.path(OUT_DIR, paste0("tracks_made_", TISSUE, ".csv")))

n_fig <- 3L + if (is.null(TP_OBJ)) 0L else 4L + if (is.null(TP_ALL_OBJ)) 0L else 2L
say("")
say(sprintf("  %d of %d genes drawn, %d figures each, PDF and PNG at %d dpi",
            sum(made$ok), nrow(made), n_fig, FIG_DPI))
say(sprintf("  output: %s", OUT_DIR))
say(sprintf("  elapsed: %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
say("")
say("  Before these go to the collaborators, compare each _byArm against its")
say("  _byAnimal. If an arm difference is carried by one animal it will be")
say("  obvious there, and with four animals per arm that is a real possibility")
say("  rather than a formality.")

if (!is.null(TP_OBJ)) {
  say("")
  say("  The same check applies twice over to the timepoint panels. Read")
  say("  _byAnimalTimepoint before _byArmTimepoint, not after: the per-animal")
  say("  panel is the only one where a before-and-after change is visible")
  say("  within an animal rather than inferred from two pooled tracks built")
  say("  from different cells.")
  by_arm_final <- table(unname(ANIMAL_ARM[paired_animals]))
  say("")
  say(sprintf("  paired animals by arm: %s",
              paste(sprintf("%s n=%d", names(by_arm_final),
                            as.integer(by_arm_final)), collapse = ", ")))
  if (any(by_arm_final < 2)) {
    say("  An arm at n=1 cannot support a treatment claim from this figure.")
    say("  It can support a description of that animal, which is a different")
    say("  sentence and worth writing as one.")
  }
} else {
  say("")
  say(sprintf("  No timepoint panels: %s has a single timepoint.", TISSUE))
}

say("")
say("  The _window_peaks.csv files say which peaks are promoter and which are")
say("  distal, so Aim 1 can be stated from a table rather than read off an image.")

invisible(gc(FALSE))
