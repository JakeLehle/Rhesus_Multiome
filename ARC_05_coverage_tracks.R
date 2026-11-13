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

# -----------------------------------------------------------------------------
# WHICH VENDOR DA CALLS COUNT AS AN OVERLAP
#
# peaks_da.csv holds 306,854 rows spanning both tissues, all three colon
# contrasts and every cell type. An unfiltered overlap against that union says
# only "this peak was called differential somewhere in the delivery", which is
# why an earlier version reported RPL18 and HSP90AB1 at 7 of 7: a union that
# large hits almost everything.
#
# So the overlap is restricted to the tissue, the cell type being plotted and
# the contrast that matches the arm comparison in the figure.
#
# These calls remain the vendor's cell-level FindMarkers results, which treat
# thousands of cells as independent replicates when the experimental unit is 8
# animals. The column is an annotation of where they looked, never evidence.
# The replicated statistic is in da_pseudobulk_<tissue>.csv from ARC_06.
# -----------------------------------------------------------------------------
VENDOR_CONTRAST <- list(Brain = "THCCBD_vs_VEH", Colon = "VEH_vs_THCCBD")

# Cell type is matched as a regex against celltype_display, not compared for
# equality, because brain microglia are three clusters (Microglia-1 at 1,489
# cells, Microglia-2 at 124, Microglia-3 at 70). Pooling them gives roughly 840
# cells per arm; splitting leaves Microglia-3 at about 35 per arm, which is too
# thin to read as a coverage track. Per-cluster counts are printed so you can
# see what went in.
PANEL <- list(
  Brain = list(
    # Microglia-1 ONLY, not ^Microglia, changed 2026-11-06.
    #
    # Two reasons. First, ARC_06 showed that 6 of 18 brain cell types have no
    # vendor RNA statistics at all and Microglia-2 and Microglia-3 are among
    # them, so pooling ^Microglia put 194 untested cells into a track built to
    # illustrate an RNA result. Second, every gene that survived selection did
    # so in Microglia-1: 14 of 20 genes have Microglia-1 as their strongest
    # cell type. Pooling three clusters to show a Microglia-1 effect dilutes it.
    celltype_regex = "^Microglia-1$",
    celltype_label = "Microglia-1",
    # The token this cell type carries in the vendor's tidy tables, which is
    # NOT the display label. Needed to filter peaks_da.csv down to the cell
    # type actually being plotted; see VENDOR_CONTRAST below.
    vendor_token = "Microglia_1",

    # The full report panel minus the two mitochondrial genes. MT-ND1 and
    # MT-ND5 are excluded because they cannot have ATAC peaks in this build at
    # all: the reference set non_nuclear_contigs = ["chrM"], cellranger-arc
    # excludes chrM from peak calling, and the consensus set came back on 22
    # contigs. Drawing them would produce 2 empty tracks, not 2 null results.
    #
    # AIM2 was in the earlier panel from the collaborators' first ask and is
    # kept, though it is not in the report's table.
    #
    # Statistics for these genes are in
    #   .../02_arc/06_selection/Brain/deg_hits_Brain.csv
    # and only 4 of them clear |log2FC| > 1 on a panel-local FDR, which is a
    # prioritisation and not an error rate. These figures are therefore
    # exploratory for most of the panel. Annotate accordingly.
    genes = c(
      # immune regulation
      "NRP1", "VSIG4", "FOLR2", "SIGLEC9", "SIGLEC5",
      # mitochondrial OXPHOS, nuclear-encoded only
      "COX4I1", "NDUFB4", "ATP5MC3",
      # lysosomal
      "ATP6V0A1", "LAMTOR2", "RNASET2",
      # ribosomal and proteostasis
      "RPL18", "RPL29", "RPS6", "RPS9", "EEF1A1", "UBA52", "NEDD8",
      "PSMB6", "HSP90AB1",
      # from the collaborators' first request
      "AIM2"
    )
  ),
  Colon = list(
    celltype_regex = "^Macrophage",
    celltype_label = "Macrophages",
    vendor_token = "Macrophage",
    # ISG15 added 2026-11-05 at the collaborators' request. It is the gut
    # counterpart of the brain report's claim that THC:CBD microglia show
    # decreased interferon-stimulated gene expression, so it tests whether that
    # ISG attenuation also appears in colon myeloid cells.
    #
    # 30 genes added 2026-11-11 from the collaboration team's colon macrophage
    # list, grouped by the biology they were presumably picked for. Three of them
    # were edited and the reason is recorded here rather than in a side note,
    # because a silent symbol fix is indistinguishable from a typo later:
    #
    #   RIG-I     -> DDX58. RIG-I is the protein name. No annotation carries it
    #                as a gene symbol, so the request as written would have
    #                missed and been reported as a resolution failure.
    #   CEMIP     -> kept, with KIAA1199 as the first alias. CEMIP is the post-
    #                2016 symbol and Mmul8 release 97 predates the rename in
    #                most of its gene models.
    #   MMP1      -> kept, with MMP1A and MMP1L1 as aliases. The macaque MMP1
    #                locus is annotated inconsistently across builds.
    #
    # Two more to watch in the resolution report rather than pre-emptively fix:
    #   FCAR  (CD89) is a pseudogene or absent in several Old World monkeys, so
    #         a miss here is biology and not an annotation gap.
    #   IGHM  sits in the immunoglobulin locus, which Mmul8 annotates poorly,
    #         and it is a B cell gene. If it resolves in macrophages, read the
    #         per-cell heatmap before believing the track: ambient Ig mRNA and
    #         doublets both land there.
    genes = c(
      # originally requested
      "HMOX1", "ISG15",
      # scavenger and efferocytosis receptors
      "MERTK", "CD36", "COLEC12", "MSR1", "MRC2", "CD14",
      # tissue remodelling and fibrosis
      "CEMIP", "PDGFRB", "TNC", "MMP1", "COL5A3", "PLAU", "F13A1",
      # interferon and antiviral
      "STAT1", "CXCL9", "DDX58", "MX1", "OAS2",
      # inflammasome, autophagy, apoptosis
      "NLRP3", "CASP9", "ULK1", "LACC1", "SERPINB1",
      # immunoregulatory and metabolic
      "IL10RA", "HSD11B1", "HIF1A",
      # adhesion and transcriptional
      "ITGA6", "PKNOX2",
      # the two to watch, see the note above
      "FCAR", "IGHM"
    )
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
  HMOX1   = c("HMOX1", "HO-1", "HMOX-1"),
  # ISG15's older human aliases are worth trying because non-human primate
  # builds sometimes carry the pre-2004 naming.
  ISG15   = c("ISG15", "G1P2", "UCRP", "ISG-15", "IFI15", "IGI15"),
  # Added with the full brain panel. ATP5MC3 was ATP5G3 before the 2016
  # nomenclature change and Mmul8 predates that, so both are tried.
  ATP5MC3 = c("ATP5MC3", "ATP5G3"),
  COX4I1  = c("COX4I1", "COX4"),
  NDUFB4  = c("NDUFB4"),
  ATP6V0A1 = c("ATP6V0A1", "ATP6V0A1", "ATP6N1", "VPP1"),
  LAMTOR2 = c("LAMTOR2", "ROBLD3", "MAPBPIP"),
  RNASET2 = c("RNASET2"),
  HSP90AB1 = c("HSP90AB1", "HSP90B", "HSPCB"),
  UBA52   = c("UBA52", "UBCEP2"),
  NEDD8   = c("NEDD8"),
  PSMB6   = c("PSMB6"),
  VSIG4   = c("VSIG4", "CRIg", "Z39IG"),
  FOLR2   = c("FOLR2", "FBP"),

  # --- colon macrophage panel, added 2026-11-11 -----------------------------
  # Older symbols first where Mmul8 release 97 predates a rename, since the
  # resolver tries these in order and the first hit wins.
  DDX58    = c("DDX58", "RIGI", "RIG-I", "RIGI1", "RLR1"),
  CEMIP    = c("KIAA1199", "CEMIP", "TMEM2L"),
  MMP1     = c("MMP1", "MMP1A", "MMP1L1", "CLG"),
  MERTK    = c("MERTK", "MER", "EYK", "NYK"),
  CD36     = c("CD36", "GP4", "GP3B", "SCARB3", "FAT"),
  COLEC12  = c("COLEC12", "CL-P1", "SCARA4", "NSR2"),
  MSR1     = c("MSR1", "CD204", "SCARA1", "MSR"),
  MRC2     = c("MRC2", "ENDO180", "CD280", "KIAA0709", "UPARAP"),
  CD14     = c("CD14"),
  PDGFRB   = c("PDGFRB", "PDGFR", "PDGFR1", "CD140B"),
  TNC      = c("TNC", "HXB", "TN"),
  COL5A3   = c("COL5A3"),
  PLAU     = c("PLAU", "UPA", "URK"),
  F13A1    = c("F13A1", "F13A"),
  STAT1    = c("STAT1"),
  CXCL9    = c("CXCL9", "MIG", "SCYB9", "CMK"),
  MX1      = c("MX1", "MXA", "IFI78"),
  OAS2     = c("OAS2"),
  NLRP3    = c("NLRP3", "NALP3", "CIAS1", "PYPAF1", "CLR1.1"),
  CASP9    = c("CASP9", "MCH6", "APAF3", "ICE-LAP6"),
  ULK1     = c("ULK1", "ATG1", "ATG1A", "UNC51"),
  # LACC1 was C13orf31 before 2014 and Mmul8 may still carry the locus tag.
  LACC1    = c("LACC1", "C13orf31", "FAMIN"),
  SERPINB1 = c("SERPINB1", "ELANH2", "MNEI", "PI2"),
  IL10RA   = c("IL10RA", "IL10R", "CD210", "CDW210A"),
  HSD11B1  = c("HSD11B1", "HSD11", "HSD11L", "HSD11B"),
  HIF1A    = c("HIF1A", "HIF1", "MOP1", "PASD8", "BHLHE78"),
  ITGA6    = c("ITGA6", "CD49F", "VLA6"),
  PKNOX2   = c("PKNOX2"),
  FCAR     = c("FCAR", "CD89", "FCAR1"),
  IGHM     = c("IGHM", "IGHM1", "MU", "VH")
)
# When every alias misses, candidates are offered from this prefix.
GENE_FAMILY_HINT <- list(
  SIGLEC9 = "SIGLEC", SIGLEC5 = "SIGLEC",
  NRP1 = "NRP", AIM2 = "AIM", HMOX1 = "HMOX", ISG15 = "ISG",
  ATP5MC3 = "ATP5", COX4I1 = "COX4", NDUFB4 = "NDUF",
  ATP6V0A1 = "ATP6V", LAMTOR2 = "LAMTOR", RNASET2 = "RNASET",
  RPL18 = "RPL", RPL29 = "RPL", RPS6 = "RPS", RPS9 = "RPS",
  EEF1A1 = "EEF1", UBA52 = "UBA", NEDD8 = "NEDD", PSMB6 = "PSMB",
  HSP90AB1 = "HSP90", VSIG4 = "VSIG", FOLR2 = "FOLR",
  # colon macrophage panel
  DDX58 = "DDX5", CEMIP = "KIAA11", MMP1 = "MMP", MERTK = "MERT",
  CD36 = "CD3", COLEC12 = "COLEC", MSR1 = "MSR", MRC2 = "MRC",
  CD14 = "CD1", PDGFRB = "PDGFR", TNC = "TN", COL5A3 = "COL5",
  PLAU = "PLAU", F13A1 = "F13", STAT1 = "STAT", CXCL9 = "CXCL",
  MX1 = "MX", OAS2 = "OAS", NLRP3 = "NLRP", CASP9 = "CASP",
  ULK1 = "ULK", LACC1 = "LACC", SERPINB1 = "SERPINB", IL10RA = "IL10R",
  HSD11B1 = "HSD11", HIF1A = "HIF1", ITGA6 = "ITGA", PKNOX2 = "PKNOX",
  FCAR = "FC", IGHM = "IGH"
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
# 0 by default: draw every track and let its n speak. Dropping a track is what
# created the NA that crashed the run, and it is also the behavior you asked me
# to stop: a control track must not disappear because it is small. Raise this
# only if you want tiny tracks gone, and the per-animal panel then subsets them
# out rather than leaving NAs in a grouping variable.
MIN_CELLS_DROP   <- 0

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

# -----------------------------------------------------------------------------
# PER-GENE WINDOW OVERRIDE, FOR GENES TOO SMALL TO READ IN THE COMMON WINDOW
#
# The 20 kb default is set by the vendor's own gene-to-peak map, which used
# -20 kb to +2 kb, so a peak outside it has no counterpart in their DA tables.
# That is the right default and it keeps genes comparable to each other.
#
# It does not suit a very small gene. ISG15 is roughly 1 to 2 kb of genomic
# span, so in a 25 kb window the gene body is a sliver at one edge and the
# promoter signal the collaborators are asking about is a few pixels wide.
# 10 kb up and 2 kb down keeps the upstream regulatory region and still leaves
# the gene legible.
#
# This is a deliberate deviation from the common window, so the figure subtitle
# says which window it used. Delete the ISG15 line to put it back on the
# default and regain strict comparability with HMOX1.
# -----------------------------------------------------------------------------
GENE_WINDOW <- list(
  ISG15 = c(up = 10000, down = 2000)
)

window_for <- function(g) {
  w <- GENE_WINDOW[[g]]
  if (is.null(w)) return(c(up = EXTEND_UP, down = EXTEND_DOWN))
  c(up = as.numeric(w[["up"]]), down = as.numeric(w[["down"]]))
}
PROMOTER_FLANK <- 2000

# -----------------------------------------------------------------------------
# LARGE GENES NEED A THIRD SCALE
#
# The byArm window is the gene body plus flanks, so for a long gene a promoter
# peak becomes invisible. NRP1 is 158,677 bp, which makes its byArm window about
# 184 kb, and chr9:33,563,188-33,564,105, the only accessibility peak in the
# whole brain dataset that reaches FDR < 0.05 on the per-animal test, is 918 bp
# of that: half a percent of the plot width.
#
# The promoter panel at TSS +/- 2 kb shows that peak but shows nothing around
# it. So genes over LARGE_GENE_BP get one more panel at TSS +/- TSS_WINDOW,
# which is wide enough to carry neighbouring regulatory elements and narrow
# enough to read. Eight of the 21 brain genes qualify.
# -----------------------------------------------------------------------------
LARGE_GENE_BP <- 20000
TSS_WINDOW    <- 10000

# -----------------------------------------------------------------------------
# RAW SIGNAL AND THE PER-CELL DISTRIBUTION
#
# WHY THIS EXISTS, which is the whole point of this revision.
#
# CoveragePlot's y axis is NOT raw Tn5 insertions. Signac divides each group's
# coverage by a per-group scale factor equal to
#
#     (number of cells in the group) x (mean sequencing depth of the group)
#
# and then multiplies by a common factor, by default the median of those group
# scale factors. So the y axis is closer to "insertions per cell per unit
# depth" than to a read count.
#
# That normalization is correct for asking "is this locus more accessible per
# cell in one arm than the other", and it is BLIND BY CONSTRUCTION to the
# question the collaborators raised. A group where 10% of cells carry all the
# signal and a group where 100% of cells carry a tenth each produce the SAME
# normalized track. Dividing by cell number erases exactly the subpopulation
# structure we want to see.
#
# So this section computes, from the fragment files, a bin-by-cell count matrix
# over the plotted window. One computation then feeds three things the
# normalized track cannot give:
#
#   1. per-group statistics: median, max, detection rate and how concentrated
#      the signal is in the top decile of cells
#   2. a RAW coverage track, summed fragments per bin, both arms on one shared
#      axis, so magnitude differences are visible rather than normalized away
#   3. a per-cell heatmap, cells sorted by total signal within arm, which shows
#      a subpopulation directly if one exists
#
# RAW_BIN_BP trades resolution against file size. 200 bp is finer than the
# 100 bp smoothing window CoveragePlot applies by default.
# -----------------------------------------------------------------------------
RAW_BIN_BP <- 200

# Track labels gain the window median and max per cell. Set FALSE if the labels
# get too long for a figure you are putting in a slide.
LABEL_WITH_CELL_STATS <- TRUE

# Cells drawn in the heatmap per group. All cells are used for the statistics
# regardless; this only caps what is rendered.
HEATMAP_MAX_CELLS <- 600

# Shared y axis across every Signac panel for one gene. The normalized value
# depends on the group scale factors, which differ between the arm panel and
# the per-animal panel, so without this the two are not on a comparable scale
# even for the same gene and the same window.
SHARE_YMAX <- TRUE

RAW_COLORS <- ARM_COLORS

FONT_SIZE <- 28          # project convention is 28 to 34
FIG_DPI   <- 300
FIG_W     <- 14
FIG_H_MIN <- 9

# -----------------------------------------------------------------------------
# THE CORNER TAG, replacing the paragraph that used to sit above every panel
#
# Changed 2026-11-12. Every figure previously carried a multi-line subtitle
# explaining the normalization, the grouping and the per-arm distribution. That
# was written to be read on its own, which is the wrong target: these figures go
# into a slide or a panel next to several others, the top line ran off the plot
# on the per-animal versions, and none of it was information a reader needs
# while looking at the picture. All of it is still in the log and in the
# _cellStats.csv next to each figure, which is where a number belongs when
# someone wants to quote it.
#
# What stays on the figure is the shortest thing that keeps it from being
# ambiguous: which tissue and cell type, and for the raw tracks the three
# numbers that change how the track should be read.
TAG_FONT_SIZE  <- FONT_SIZE      # same size as the rest, not a footnote
TAG_COLOR      <- "#333333"

# The concentration fraction, used BOTH by cell_stats() to compute the share and
# by the figure text that reports it. One constant rather than a hardcoded 0.10
# in the statistic and a hardcoded "10%" in the label, because those two drift
# apart silently and the figure would then state a number it did not measure.
TOP_FRACTION   <- 0.10

# Fill under the raw coverage traces. Two semi-transparent fills on one axis
# make the exclusive and shared parts of the signal separable at a glance.
#
# Note on what the overlap looks like: VEH is #0072B2 and THC/CBD is #D55E00,
# which are blue and vermillion, so where they overlap they blend to a dark
# olive-brown rather than purple. Purple comes from blue plus magenta. If the
# blend matters more than keeping the current arm colors, set
# RAW_COLORS <- c("VEH" = "#0072B2", "THC/CBD" = "#CC79A7") below, which is
# still Okabe-Ito and still colorblind-safe, and the overlap reads violet.
RAW_FILL_ALPHA <- 0.22

# Signac's default. Scales each group's coverage by its total reads in the
# assay, which is what makes the per-animal panel necessary rather than
# optional: it corrects depth between groups, not between animals pooled
# inside one group.
NORMALIZATION <- "RC"

# Signac 1.17.1's CoveragePlot does not accept a normalization argument, so the
# value above is dropped by call_supported() and Signac's own default applies.
# Claiming "RC normalization" in a figure subtitle when the argument was never
# passed is a false statement on the figure, so the label is resolved from the
# function's actual formals rather than from the constant.
NORM_LABEL <- if ("normalization" %in% names(formals(Signac::CoveragePlot))) {
  paste0(NORMALIZATION, " normalization")
} else {
  "Signac default normalization"
}

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
  # trim = TRUE matters. format() on a numeric VECTOR pads every element to a
  # common width, so format(c(56, 357)) gives " 56" and "357". That space then
  # ends up inside the label as "(n= 56)", which broke both strip_n() below and
  # every palette lookup built on it, silently.
  paste0(base, " (n=", format(n, big.mark = ",", trim = TRUE), ")")
}
# Tolerant of stray whitespace either way, so a padded label can never again
# slip past this and take the color mapping down with it.
strip_n <- function(x) base::sub("[[:space:]]*\\(n=[[:space:]]*[0-9,]+\\)[[:space:]]*$", "", x)

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
if (length(GENE_WINDOW)) {
  for (gw in names(GENE_WINDOW)) {
    say(sprintf("           %s overrides this: -%s to +%s bp", gw,
                format(GENE_WINDOW[[gw]][["up"]], big.mark = ","),
                format(GENE_WINDOW[[gw]][["down"]], big.mark = ",")))
  }
}
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
TP_AN_OBJ <- NULL
TP_ALL_OBJ <- NULL
paired_animals <- character(0)
if (!HAS_TIMECOURSE) {
  subh("timepoint panels")
  # tp_levels is the intersection with TIMEPOINT_ORDER, which excludes
  # "Single timepoint", so for brain it is empty and printed as "()".
  tp_shown <- if (length(tp_levels)) tp_levels else levels(sub_obj$tp)
  wrn(paste0(TISSUE, " has a single timepoint (", paste(tp_shown, collapse = ", "),
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
        n <- sum(as.character(o$arm_tp) == nm, na.rm = TRUE)
        say(sprintf("  %-40s %6d cells%s", nm, n,
                    if (n < MIN_CELLS_FLAG) "   (thin)" else ""))
      }
      # The object carries arm as its own column, so read it. Reconstructing it
      # from "VEH Pre-infection (n= 56)" by regex is what produced a false
      # NO CONTROL TRACK warning on a panel that plainly had a control.
      present_arms <- unique(as.character(o$arm[!is.na(o$arm)]))
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
        n <- sum(as.character(o$animal_tp) == nm, na.rm = TRUE)
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

    # The per-animal panel gets its own object with the dropped cells removed,
    # rather than being handed a grouping variable containing NAs. With
    # MIN_CELLS_DROP at 0 nothing is dropped and this is the same object; it
    # matters only if you raise that threshold.
    an_keep <- colnames(TP_OBJ)[!is.na(TP_OBJ$animal_tp)]
    TP_AN_OBJ <- if (length(an_keep) == ncol(TP_OBJ)) TP_OBJ else {
      wrn(sprintf("per-animal panel excludes %d cell(s) in dropped tracks",
                  ncol(TP_OBJ) - length(an_keep)))
      o <- subset(TP_OBJ, cells = an_keep)
      o$animal_tp <- droplevels(o$animal_tp)
      o
    }
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
  # peaks_da.csv from MULTIOME_01 names its peak column peak_id and already
  # carries chrom, start and end as separate columns. The earlier version of
  # this block looked for atac_peak / peak / region and then regex-parsed a
  # coordinate string, so it failed its tryCatch on every run and the vendor DA
  # annotation was silently absent from every _window_peaks.csv ever written.
  # Verified against MULTIOME_01_build_tidy.py cell 4.
  da_ok <- tryCatch({
    da <- fread(DA_PEAKS_CSV)
    if (all(c("chrom", "start", "end") %in% names(da))) {
      good <- !is.na(da$chrom) & !is.na(da$start) & !is.na(da$end)
      gr <- GRanges(as.character(da$chrom[good]),
                    IRanges(as.numeric(da$start[good]), as.numeric(da$end[good])))
    } else {
      # Fallback for a table that only has the composite id.
      pk_col <- intersect(c("peak_id", "atac_peak", "peak", "region"), names(da))[1]
      if (is.na(pk_col)) {
        stop("no chrom/start/end and no peak id column; columns are: ",
             paste(names(da), collapse = ", "))
      }
      pv <- as.character(da[[pk_col]])
      m <- regmatches(pv, regexec("^(.*)-(\\d+)-(\\d+)$", pv))
      good <- vapply(m, length, integer(1)) == 4
      gr <- GRanges(vapply(m[good], `[`, character(1), 2),
                    IRanges(as.numeric(vapply(m[good], `[`, character(1), 3)),
                            as.numeric(vapply(m[good], `[`, character(1), 4))))
    }
    ct_col <- intersect(c("celltype", "atac_celltype", "cell_type"), names(da))[1]
    if (!is.na(ct_col)) gr$celltype <- as.character(da[[ct_col]])[good]
    ts_col <- intersect(c("tissue"), names(da))[1]
    if (!is.na(ts_col)) gr$tissue <- as.character(da[[ts_col]])[good]
    cn_col <- intersect(c("contrast"), names(da))[1]
    if (!is.na(cn_col)) gr$contrast <- as.character(da[[cn_col]])[good]
    gr
  }, error = function(e) { wrn(paste0("could not use ", DA_PEAKS_CSV, ": ", conditionMessage(e))); NULL })
  if (!is.null(da_ok)) {
    n_all <- length(da_ok)
    tok <- PANEL[[TISSUE]]$vendor_token
    con <- VENDOR_CONTRAST[[TISSUE]]
    keep <- rep(TRUE, n_all)
    if (!is.null(da_ok$tissue))   keep <- keep & da_ok$tissue == TISSUE
    if (!is.null(da_ok$celltype) && !is.null(tok)) keep <- keep & da_ok$celltype == tok
    if (!is.null(da_ok$contrast) && !is.null(con)) keep <- keep & da_ok$contrast == con
    keep[is.na(keep)] <- FALSE
    da_gr <- da_ok[keep]
    say(sprintf("  vendor DA peaks: %s loaded, %s after filtering to %s / %s / %s",
                format(n_all, big.mark = ","), format(length(da_gr), big.mark = ","),
                TISSUE, if (is.null(tok)) "any cell type" else tok,
                if (is.null(con)) "any contrast" else con))
    if (!length(da_gr)) {
      wrn("no vendor DA peaks survive that filter; the overlap column will be")
      say("          all FALSE. Check the token and contrast against")
      say("          peaks_da.csv before reading anything into it.")
      da_gr <- NULL
    } else {
      ok("overlap column reflects THIS cell type and contrast only")
      say("          (vendor calls are cell-level and pseudoreplicated: an")
      say("          annotation of where they looked, not evidence)")
    }
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

# Counts what THIS run wrote. The previous version counted every PDF in the
# output directory, so a run whose panels all failed still reported 121 files,
# most of them stale from an earlier run.
FIGS_WRITTEN <- new.env(parent = emptyenv())
FIGS_WRITTEN$n <- 0L
FIGS_WRITTEN$stems <- character(0)

save_fig <- function(p, stem, n_tracks) {
  h <- max(FIG_H_MIN, 3 + 1.6 * n_tracks)
  pdf_path <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(stem, ".png"))
  ggsave(pdf_path, p, width = FIG_W, height = h, units = "in",
         device = grDevices::cairo_pdf, limitsize = FALSE)
  ggsave(png_path, p, width = FIG_W, height = h, units = "in",
         dpi = FIG_DPI, limitsize = FALSE)
  FIGS_WRITTEN$n <- FIGS_WRITTEN$n + 1L
  FIGS_WRITTEN$stems <- c(FIGS_WRITTEN$stems, stem)
  say(sprintf("    wrote %s.pdf and .png  (%.0f x %.0f in, %d dpi)",
              stem, FIG_W, h, FIG_DPI))
}

# The tag text. A function rather than a constant so it reads TISSUE and
# CT_LABEL at call time, after Phase 2 has set them.
#
#   tag()              -> "Colon Macrophages"
#   tag("by animal")   -> "Colon Macrophages  |  by animal"
#
# Nothing else goes in. The grouping token is there only because byArm, byAnimal
# and the promoter zoom are otherwise indistinguishable once the figure is out
# of its filename, and it is two or three words in every case.
tag <- function(what = NULL) {
  base <- paste(TISSUE, CT_LABEL)
  if (is.null(what) || !nzchar(what)) base else paste0(base, "  |  ", what)
}

# Gene names are italic, everything else roman. The gene goes in the title on
# its own so the whole title can be italicized. The tag is right-aligned so it
# sits in the top right of the figure and cannot collide with the gene name or
# run past the edge, which is what the old multi-line subtitle did on the
# per-animal panels.
decorate <- function(p, gene, tag_text) {
  p +
    plot_annotation(
      title = gene,
      subtitle = tag_text,
      theme = theme(
        plot.title    = element_text(size = FONT_SIZE + 4, face = "bold.italic"),
        plot.subtitle = element_text(size = TAG_FONT_SIZE, hjust = 1,
                                     colour = TAG_COLOR)
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

draw <- function(object, region, group_by, cols, extend_up, extend_down,
                 ymax = NULL) {
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
      ymax              = ymax,
      normalization     = NORMALIZATION
    ),
    label = "CoveragePlot")
  if (!"cols" %in% names(formals(CoveragePlot))) p <- apply_arm_colors(p, cols)
  p
}

# -----------------------------------------------------------------------------
# Raw signal, per-cell statistics, and the subpopulation view
# -----------------------------------------------------------------------------

# The plotted window, strand-aware, matching what draw() shows.
plotted_window <- function(region, extend_up, extend_down) {
  w <- region
  if (as.character(strand(w)) == "-") {
    start(w) <- max(1, start(region) - extend_down); end(w) <- end(region) + extend_up
  } else {
    start(w) <- max(1, start(region) - extend_up);   end(w) <- end(region) + extend_down
  }
  w
}

window_bins <- function(w, bin = RAW_BIN_BP) {
  st <- seq(from = start(w), to = end(w), by = bin)
  GRanges(as.character(seqnames(w))[1],
          IRanges(st, pmin(st + bin - 1L, end(w))))
}

# Bins x cells of raw fragment counts, straight off the fragment files. This is
# the one computation behind the raw track, the statistics and the heatmap.
window_counts <- function(object, w) {
  bins <- window_bins(w)
  if (!length(bins)) return(NULL)
  fm <- tryCatch(
    call_supported(
      FeatureMatrix,
      required = list(fragments = Fragments(object[["ATAC"]]),
                      features  = bins,
                      cells     = colnames(object)),
      optional = list(process_n = 20000, verbose = FALSE),
      label = "FeatureMatrix"),
    error = function(e) { wrn(paste0("window_counts failed: ", conditionMessage(e))); NULL })
  fm
}

# Per-group distribution statistics. Every number here is a RAW count, so none
# of it is affected by the normalization discussed at the top of this file.
#
# top_decile_share is the diagnostic for the collaborators' question: the
# fraction of all window fragments contributed by the most-covered 10% of
# cells. At 0.1 the signal is spread evenly; near 1.0 it is a handful of cells.
cell_stats <- function(fm, grp) {
  if (is.null(fm)) return(NULL)
  per_cell <- Matrix::colSums(fm)
  grp <- as.character(grp[colnames(fm)])
  out <- rbindlist(lapply(sort(unique(grp[!is.na(grp)])), function(gg) {
    v <- per_cell[which(grp == gg)]
    if (!length(v)) return(NULL)
    k <- max(1L, ceiling(TOP_FRACTION * length(v)))
    top <- sort(v, decreasing = TRUE)[seq_len(k)]
    data.table(
      group            = gg,
      n_cells          = length(v),
      total_fragments  = sum(v),
      mean_per_cell    = mean(v),
      median_per_cell  = stats::median(v),
      q90_per_cell     = unname(stats::quantile(v, 0.90)),
      max_per_cell     = max(v),
      pct_cells_with_any = 100 * mean(v > 0),
      top_decile_share = if (sum(v) > 0) sum(top) / sum(v) else NA_real_,
      frag_per_cell_per_kb = mean(v) / (nrow(fm) * RAW_BIN_BP / 1000)
    )
  }))
  out
}

# Label gains the raw median and max, which is what distinguishes "every cell a
# little" from "a few cells a lot".
label_with_stats <- function(base_levels, st) {
  if (!LABEL_WITH_CELL_STATS || is.null(st)) return(setNames(base_levels, base_levels))
  out <- vapply(base_levels, function(lv) {
    r <- st[group == lv]
    if (!nrow(r)) return(lv)
    sprintf("%s  med %g, max %g", lv, r$median_per_cell[1], r$max_per_cell[1])
  }, character(1))
  setNames(out, base_levels)
}

# The raw track. Both arms on ONE shared axis in a single panel, because the
# question is a magnitude comparison and separate panels invite separate scales.
raw_track_plot <- function(fm, grp, w, gene, tag_text, st) {
  if (is.null(fm)) return(NULL)
  bins <- window_bins(w)
  gv <- as.character(grp[colnames(fm)])
  keep <- !is.na(gv)
  lv <- intersect(names(ARM_COLORS), unique(gv[keep]))
  if (!length(lv)) lv <- sort(unique(gv[keep]))
  # TWO VERTICES PER BIN, AT THE BIN EDGES, NOT ONE AT THE MIDPOINT
  #
  # Changed 2026-11-13. The previous version put one vertex at each bin's
  # midpoint, then drew geom_step for the outline and geom_area for the fill.
  # Those two geoms disagree about what happens between vertices: geom_step
  # holds the value flat and then jumps, while geom_area interpolates in a
  # straight diagonal. So on every rising or falling edge the fill cut a
  # triangle off under the line and stopped short of it, which is exactly the
  # "shaded region doesn't push all the way up to the lines" that showed up in
  # the figures. It was worst on sharp peaks, where the error is largest.
  #
  # Emitting a vertex at the start and the end of each bin makes the data
  # itself the step: a plain geom_line through those vertices IS the staircase,
  # and a geom_ribbon under the same vertices fills it exactly. The fill and the
  # outline are now the same shape by construction rather than by two geoms
  # agreeing, so they cannot drift apart again.
  #
  # It also puts the trace where the data actually is. A bin covering
  # 1,000-1,200 bp is now drawn across 1,000-1,200 bp instead of being centred
  # on 1,100 and spilling half a bin past each end of the window.
  edges <- as.vector(rbind(start(bins) - 1L, end(bins)))
  dt <- rbindlist(lapply(lv, function(gg) {
    cols <- which(gv == gg)
    cnt <- as.numeric(Matrix::rowSums(fm[, cols, drop = FALSE]))
    data.table(group = gg, pos = edges, count = rep(cnt, each = 2L))
  }))
  dt[, group := factor(group, levels = lv)]

  # One line per arm, carrying the three numbers that change how the trace
  # should be read: n says how much cell is behind it, the read median says
  # whether a
  # typical cell contributes at all, and the decile share says whether the trace
  # is an average over cells that mostly contribute nothing.
  #
  # Tissue and cell type go in the right-aligned subtitle rather than in here.
  # Four lines of 28 pt inside the panel is tall enough to land on a peak, and
  # two is not.
  #
  # The decile is fixed at 10% by cell_stats(), so the first number in the
  # phrase is a constant rather than something read off the data. It is written
  # out anyway because "top 10% = 42%" does not say which 10% or a percentage
  # of what, and the figure goes to readers who did not build it.
  cap <- st[match(lv, st$group)]
  cap <- cap[!is.na(cap$group)]
  corner <- paste(sprintf(
    "%s  n=%d  read median %g  top %.0f%% of cells hold %.0f%% of reads",
    cap$group, cap$n_cells, cap$median_per_cell,
    100 * TOP_FRACTION, 100 * cap$top_decile_share), collapse = "\n")

  ggplot(dt, aes(x = pos, y = count)) +
    # Fill first, outline over it, both from the same bin-edge vertices.
    #
    # geom_ribbon rather than geom_area because geom_area stacks by default and
    # would make the second arm sit on top of the first, so the figure would
    # read as a sum rather than a comparison. A ribbon from 0 is the same shape
    # with no stacking to switch off.
    geom_ribbon(aes(ymin = 0, ymax = count, fill = group),
                alpha = RAW_FILL_ALPHA, colour = NA) +
    # geom_line, not geom_step. The staircase is already in the vertices, so
    # geom_step here would step a second time and offset the outline by half a
    # bin from its own fill.
    geom_line(aes(colour = group), linewidth = 1.1) +
    scale_fill_manual(values = RAW_COLORS[lv], name = NULL, drop = FALSE) +
    scale_colour_manual(values = RAW_COLORS[lv], name = NULL, drop = FALSE) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", trim = TRUE)) +
    # Flush to the inside of the top-right panel corner. hjust and vjust are
    # exactly 1 rather than a fraction past it, because anything over 1 pushes
    # the text past the panel edge where clipping cuts it off.
    annotate("text", x = Inf, y = Inf, label = corner,
             hjust = 1, vjust = 1, size = TAG_FONT_SIZE / .pt,
             colour = TAG_COLOR, lineheight = 1.1) +
    # Headroom so the two corner lines sit above the tallest peak rather than
    # on it. 0.28 is enough for two lines at 28 pt in a 9 inch figure.
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    # One legend, not two. Fill and colour carry the same information and
    # whether ggplot merges their guides depends on version, so the fill guide
    # is dropped outright instead of being left to chance.
    guides(fill = "none") +
    labs(title = gene, subtitle = tag_text,
         x = paste0(as.character(seqnames(w))[1], " position (bp)"),
         y = sprintf("raw fragments per %d bp", RAW_BIN_BP)) +
    theme_minimal(base_size = FONT_SIZE) +
    theme(plot.title = element_text(face = "bold.italic", size = FONT_SIZE + 4),
          plot.subtitle = element_text(size = TAG_FONT_SIZE, hjust = 1,
                                       colour = TAG_COLOR),
          legend.position = "top",
          panel.grid.minor = element_blank())
}

# ONE ROW PER CELL. x is genomic position in RAW_BIN_BP bins, fill is the raw
# fragment count that one nucleus contributed to that bin, and rows are ordered
# by each cell's total over the whole window, highest at the bottom of the rank
# axis. Rank restarts within each arm, so row 1 of VEH and row 1 of THC/CBD are
# different cells and the two facets are not aligned row for row.
#
# Only non-zero bins are drawn, so a cell with nothing in the window renders as
# an unbroken strip of panel background. That is the point: if a subpopulation
# carries the locus it shows as a dense block at the top of one arm with empty
# background filling the rest, whereas a broad shift shows as speckle running
# the full height.
#
# HEATMAP_MAX_CELLS caps the rows per arm so the raster stays legible. The cap
# takes the HIGHEST-signal cells, so a capped facet is biased toward coverage
# and the y axis label says when that happened.
percell_heatmap <- function(fm, grp, w, gene, tag_text) {
  if (is.null(fm)) return(NULL)
  bins <- window_bins(w)
  gv <- as.character(grp[colnames(fm)])
  lv <- intersect(names(ARM_COLORS), unique(gv[!is.na(gv)]))
  if (!length(lv)) lv <- sort(unique(gv[!is.na(gv)]))
  per_cell <- Matrix::colSums(fm)
  rows <- list()
  capped <- FALSE
  n_tot <- integer(0)
  for (gg in lv) {
    cols <- which(gv == gg)
    n_tot[gg] <- length(cols)
    cols <- cols[order(per_cell[cols], decreasing = TRUE)]
    if (length(cols) > HEATMAP_MAX_CELLS) {
      cols <- cols[seq_len(HEATMAP_MAX_CELLS)]
      capped <- TRUE
    }
    sub <- fm[, cols, drop = FALSE]
    nz <- Matrix::summary(as(sub, "TsparseMatrix"))
    if (!nrow(nz)) next
    rows[[gg]] <- data.table(
      group = gg,
      pos   = ((start(bins) + end(bins)) / 2)[nz$i],
      cell_rank = nz$j,
      count = nz$x)
  }
  if (!length(rows)) return(NULL)
  dt <- rbindlist(rows)
  dt[, group := factor(group, levels = lv)]
  ylab <- if (capped) {
    sprintf("cell, ranked by signal (top %s per arm of %s)",
            format(HEATMAP_MAX_CELLS, big.mark = ",", trim = TRUE),
            paste(format(n_tot[lv], big.mark = ",", trim = TRUE), collapse = " / "))
  } else "cell, ranked by signal"
  ggplot(dt, aes(x = pos, y = cell_rank, fill = count)) +
    geom_raster() +
    facet_wrap(~ group, ncol = 1, scales = "free_y") +
    scale_fill_gradient(low = "#DCE6F1", high = "#B2182B", name = "frags",
                        trans = "sqrt") +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", trim = TRUE)) +
    labs(title = gene, subtitle = tag_text,
         x = paste0(as.character(seqnames(w))[1], " position (bp)"),
         y = ylab) +
    theme_minimal(base_size = FONT_SIZE) +
    theme(plot.title = element_text(face = "bold.italic", size = FONT_SIZE + 4),
          plot.subtitle = element_text(size = TAG_FONT_SIZE, hjust = 1,
                                       colour = TAG_COLOR),
          strip.text = element_text(size = FONT_SIZE, face = "bold"),
          # A near-white panel so "this cell has nothing here" is visibly
          # different from the pale blue of a single fragment.
          panel.background = element_rect(fill = "#FAFAFA", colour = NA),
          panel.grid = element_blank())
}

# The y limit CoveragePlot actually computed, harvested from the built plot so
# it can be reused as a common ymax across every panel for one gene.
harvest_ymax <- function(p) {
  tryCatch({
    b <- ggplot2::ggplot_build(p)
    rr <- unlist(lapply(b$layout$panel_scales_y, function(sc) sc$range$range))
    v <- suppressWarnings(max(as.numeric(rr), na.rm = TRUE))
    if (is.finite(v) && v > 0) v else NULL
  }, error = function(e) NULL)
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

  # Both of these used to be printed on every figure. The figures now carry a
  # short corner tag instead, so they are logged rather than lost: the paired
  # distinction is otherwise visible only in the _allAnimals filename suffix,
  # and the normalization label was the whole point of the 2026-11-11 fix.
  tp_note <- if (PAIRED_ONLY) "paired animals only" else "all available animals"
  if (!is.null(TP_OBJ)) say(sprintf("  timepoint panels: %s", tp_note))
  say(sprintf("  Signac panels: %s", NORM_LABEL))

  # Window for THIS gene: the common one unless GENE_WINDOW overrides it.
  ew <- window_for(g)
  win_note <- if (is.null(GENE_WINDOW[[g]])) "" else
    sprintf(", window -%s/+%s bp", format(ew[["up"]], big.mark = ","),
            format(ew[["down"]], big.mark = ","))
  if (nzchar(win_note)) {
    say(sprintf("  non-default window: -%s to +%s bp (small gene)",
                format(ew[["up"]], big.mark = ","),
                format(ew[["down"]], big.mark = ",")))
  }

  raw_done <- FALSE
  res <- tryCatch({
    # --- raw signal first, because it decides the labels ---------------------
    # Done before any Signac panel so the per-cell median and max can go into
    # the track labels, and so the raw figures exist even if CoveragePlot
    # struggles with an awkward window.
    wplot <- plotted_window(region, ew[["up"]], ew[["down"]])
    fm <- window_counts(POST_OBJ, wplot)
    arm_vec <- setNames(as.character(POST_OBJ$arm), colnames(POST_OBJ))
    st <- cell_stats(fm, arm_vec)

    if (!is.null(st)) {
      say(sprintf("  raw signal over %s bp in %d bins of %d bp",
                  format(width(wplot), big.mark = ","),
                  length(window_bins(wplot)), RAW_BIN_BP))
      for (k in seq_len(nrow(st))) {
        r <- st[k]
        say(sprintf("    %-10s n=%4d  total %6.0f  median %4g  q90 %5g  max %5g  %5.1f%% covered  top decile %4.0f%%",
                    r$group, r$n_cells, r$total_fragments, r$median_per_cell,
                    r$q90_per_cell, r$max_per_cell, r$pct_cells_with_any,
                    100 * r$top_decile_share))
      }
      conc <- st[!is.na(top_decile_share) & top_decile_share > 0.5]
      if (nrow(conc)) {
        wrn(sprintf("%s: the top 10%% of cells hold over half the signal, so the",
                    paste(conc$group, collapse = " and ")))
        say("          normalized track for that arm is an average over cells that")
        say("          mostly contribute nothing. Read the heatmap before the track.")
      }
      st_out <- copy(st)
      st_out[, `:=`(gene = g, tissue = TISSUE, celltype = CT_LABEL,
                    window = sprintf("%s:%d-%d", as.character(seqnames(wplot))[1],
                                     start(wplot), end(wplot)))]
      fwrite(st_out, file.path(OUT_DIR, paste0(stem_base, "_cellStats.csv")))
      say(sprintf("    wrote %s_cellStats.csv", stem_base))

      pr <- raw_track_plot(fm, arm_vec, wplot, g, tag("raw fragments"), st)
      if (!is.null(pr)) save_fig(pr, paste0(stem_base, "_rawCoverage"), 3)
      ph <- percell_heatmap(fm, arm_vec, wplot, g, tag())
      if (!is.null(ph)) save_fig(ph, paste0(stem_base, "_perCellHeatmap"), 7)
      raw_done <- TRUE
    }

    # Per-gene track labels carrying this window's raw median and max.
    #
    # arm_map is named by arm ("VEH", "THC/CBD"), and current R keeps those
    # names through factor(). Seurat's $<- then tries to match them against
    # cell barcodes, finds no overlap, and errors with "No cell overlap between
    # new meta data and Seurat object". That killed every coverage panel for
    # every gene while the raw figures above succeeded, so the failure looked
    # like a plotting problem rather than a metadata one.
    #
    # Named by barcode, which is what Seurat actually wants.
    arm_map <- label_with_stats(levels(POST_OBJ$arm), st)
    lab_vec <- unname(arm_map[as.character(POST_OBJ$arm)])
    names(lab_vec) <- colnames(POST_OBJ)
    POST_OBJ$arm_win <- factor(lab_vec,
                               levels = unname(arm_map[levels(POST_OBJ$arm)]))

    # Refuse to draw on a broken grouping column rather than emitting 21
    # identical failures. Checks what actually matters: every cell labelled,
    # and the right number of distinct tracks.
    if (!"arm_win" %in% names(POST_OBJ@meta.data) ||
        anyNA(POST_OBJ$arm_win) ||
        nlevels(POST_OBJ$arm_win) != nlevels(POST_OBJ$arm)) {
      stop("arm_win was not built correctly: ",
           sum(is.na(POST_OBJ$arm_win)), " unlabelled cells, ",
           nlevels(POST_OBJ$arm_win), " levels where ",
           nlevels(POST_OBJ$arm), " were expected")
    }
    arm_pal <- setNames(unname(ARM_COLORS[levels(POST_OBJ$arm)]),
                        unname(arm_map[levels(POST_OBJ$arm)]))
    arm_pal <- arm_pal[!is.na(arm_pal)]

    # --- the arm panels, post-infection cells only ---------------------------
    # Control track first in all three, because the figure is read as treated
    # against control and not the other way round.
    # 1. by arm, wide. Drawn first so its y limit becomes the common ymax.
    p1 <- draw(POST_OBJ, region, "arm_win", arm_pal, ew[["up"]], ew[["down"]])
    ymax_common <- if (SHARE_YMAX) harvest_ymax(p1) else NULL
    if (!is.null(ymax_common)) {
      say(sprintf("  shared ymax across this gene's panels: %.3g", ymax_common))
    }
    p1 <- decorate(p1, g, tag("by arm"))
    save_fig(p1, paste0(stem_base, "_byArm"), length(levels(POST_OBJ$arm_win)))

    # 2. by animal, same window. NOTE: Signac rescales per group, and the
    # per-animal groups have different scale factors from the arm groups, so
    # forcing the arm ymax here is what makes the two panels comparable.
    p2 <- draw(POST_OBJ, region, "animal_track", NULL, ew[["up"]], ew[["down"]],
               ymax = ymax_common)
    p2 <- decorate(p2, g, tag("by animal"))
    save_fig(p2, paste0(stem_base, "_byAnimal"), length(at))

    # 3. promoter zoom, by arm
    p3 <- draw(POST_OBJ, prom, "arm_win", arm_pal, 0, 0, ymax = ymax_common)
    p3 <- decorate(p3, g, tag(sprintf("TSS +/- %s kb",
                                      format(PROMOTER_FLANK / 1000, trim = TRUE))))
    save_fig(p3, paste0(stem_base, "_promoter"), length(levels(POST_OBJ$arm_win)))

    # 3b. TSS window, by arm, for genes too long to read at gene-body scale.
    if (as.numeric(width(region)) > LARGE_GENE_BP) {
      tssw <- region
      start(tssw) <- max(1, t0 - TSS_WINDOW); end(tssw) <- t0 + TSS_WINDOW
      strand(tssw) <- strand(region)
      say(sprintf("  gene body is %s bp, adding a TSS +/- %s bp panel",
                  format(as.numeric(width(region)), big.mark = ","),
                  format(TSS_WINDOW, big.mark = ",")))
      p3b <- draw(POST_OBJ, tssw, "arm_win", arm_pal, 0, 0, ymax = ymax_common)
      p3b <- decorate(p3b, g, tag(sprintf("TSS +/- %s kb",
                                          format(TSS_WINDOW / 1000, trim = TRUE))))
      save_fig(p3b, paste0(stem_base, "_tss", TSS_WINDOW / 1000, "kb"),
               length(levels(POST_OBJ$arm_win)))
    }

    # --- the timepoint panels, colon only ------------------------------------
    if (!is.null(TP_OBJ)) {
      # 4. pre vs post, pooled across arms. The plainest version of the
      # question the collaborators asked, and the one with the most cells
      # behind each track.
      p4 <- draw(TP_OBJ, region, "tp_lab",
                 palette_for(levels(TP_OBJ$tp_lab), TP_COLORS),
                 ew[["up"]], ew[["down"]])
      p4 <- decorate(p4, g, tag("pre vs post"))
      save_fig(p4, paste0(stem_base, "_byTimepoint"), length(levels(TP_OBJ$tp)))

      # 5. arm by timepoint. Four tracks in control-first order: VEH before,
      # VEH after, THC/CBD before, THC/CBD after. This is the panel the
      # collaborators asked for, and the one where the 3:1 split in which
      # animals have a pre-infection sample bites hardest.
      p5 <- draw(TP_OBJ, region, "arm_tp",
                 palette_for(levels(TP_OBJ$arm_tp), ARM_TP_COLORS),
                 ew[["up"]], ew[["down"]])
      p5 <- decorate(p5, g, tag("arm x timepoint"))
      save_fig(p5, paste0(stem_base, "_byArmTimepoint"),
               length(levels(TP_OBJ$arm_tp)))

      # 6. per animal, pairs adjacent. The only panel in which a within-animal
      # change is visible as such, rather than inferred from two pooled tracks
      # made of different cells.
      p6 <- draw(TP_AN_OBJ, region, "animal_tp", NULL, ew[["up"]], ew[["down"]])
      p6 <- decorate(p6, g, tag("animal x timepoint"))
      save_fig(p6, paste0(stem_base, "_byAnimalTimepoint"),
               length(levels(TP_AN_OBJ$animal_tp)))

      # 7. promoter zoom, arm by timepoint
      p7 <- draw(TP_OBJ, prom, "arm_tp",
                 palette_for(levels(TP_OBJ$arm_tp), ARM_TP_COLORS), 0, 0)
      p7 <- decorate(p7, g, tag(sprintf("TSS +/- %s kb, arm x timepoint",
                                        format(PROMOTER_FLANK / 1000, trim = TRUE))))
      save_fig(p7, paste0(stem_base, "_promoter_byArmTimepoint"),
               length(levels(TP_OBJ$arm_tp)))
    }

    # 8 and 9. The all-animals variant of the arm-by-timepoint comparison, with
    # more cells behind the control arm and the different-animal-sets confound
    # in its place. Named so it cannot be confused with the paired version.
    if (!is.null(TP_ALL_OBJ)) {
      p8 <- draw(TP_ALL_OBJ, region, "arm_tp",
                 palette_for(levels(TP_ALL_OBJ$arm_tp), ARM_TP_COLORS),
                 ew[["up"]], ew[["down"]])
      p8 <- decorate(p8, g, tag("arm x timepoint, all animals"))
      save_fig(p8, paste0(stem_base, "_byArmTimepoint_allAnimals"),
               length(levels(TP_ALL_OBJ$arm_tp)))

      p9 <- draw(TP_ALL_OBJ, prom, "arm_tp",
                 palette_for(levels(TP_ALL_OBJ$arm_tp), ARM_TP_COLORS), 0, 0)
      p9 <- decorate(p9, g, tag(sprintf("TSS +/- %s kb, all animals",
                                        format(PROMOTER_FLANK / 1000, trim = TRUE))))
      save_fig(p9, paste0(stem_base, "_promoter_byArmTimepoint_allAnimals"),
               length(levels(TP_ALL_OBJ$arm_tp)))
    }

    # peak table for the window
    wt <- window_peaks_table(symbol, region, ew[["up"]], ew[["down"]])
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

  # Two flags, not one. The raw figures and the Signac panels fail
  # independently, and the last run wrote 42 good raw figures while reporting
  # "0 of 21 genes drawn" because a single ok column collapsed both.
  made <- rbind(made, data.table(gene = g, symbol = symbol,
                                 region = paste0(as.character(seqnames(region)), ":",
                                                 start(region), "-", end(region)),
                                 raw_ok = isTRUE(raw_done),
                                 panels_ok = res,
                                 ok = res))
}

# -----------------------------------------------------------------------------
# Cell 8: summary
# -----------------------------------------------------------------------------

hdr("Done")
print(made)
fwrite(made, file.path(OUT_DIR, paste0("tracks_made_", TISSUE, ".csv")))

# Counted from what was actually written, not predicted. The TSS panel is
# conditional on gene length, so 4 of the 21 brain genes get one more figure
# than the rest and the old fixed arithmetic under-reported them.
n_fig_base <- 3L + if (is.null(TP_OBJ)) 0L else 4L + if (is.null(TP_ALL_OBJ)) 0L else 2L
n_written <- FIGS_WRITTEN$n
n_stale <- length(list.files(OUT_DIR, pattern = "\\.pdf$")) - n_written
n_fig <- n_fig_base
say("")
say(sprintf("  coverage panels: %d of %d genes   raw figures: %d of %d",
            sum(made$panels_ok), nrow(made), sum(made$raw_ok), nrow(made)))
if (sum(made$panels_ok) < nrow(made)) {
  wrn(sprintf("%d gene(s) produced no coverage panel. Any _byArm, _byAnimal,",
              nrow(made) - sum(made$panels_ok)))
  say("          _promoter or _tss file for those genes in the output directory")
  say("          is from an EARLIER run. Do not read it as this run's output.")
}
say(sprintf("  %d figures per gene (%d for genes over %s bp), plus rawCoverage",
            n_fig_base, n_fig_base + 1L, format(LARGE_GENE_BP, big.mark = ",")))
say("  and perCellHeatmap for every gene")
say(sprintf("  %d PDF files written BY THIS RUN, each with a matching PNG at %d dpi",
            n_written, FIG_DPI))
if (n_stale > 0) {
  wrn(sprintf("%d other PDF(s) in that directory are from an earlier run.", n_stale))
  say("          They are NOT regenerated figures. Check timestamps before")
  say("          sending anything from there.")
}
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
