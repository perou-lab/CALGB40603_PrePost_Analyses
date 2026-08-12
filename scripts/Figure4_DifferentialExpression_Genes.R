#!/usr/bin/env Rscript
# =====================================================================================
# Figure 4: Paired Differential Gene Expression — Pre- vs. Posttreatment,
#           Patients Who Remained Basal Subtype Throughout (Basal-Basal) (N = 21 pairs)
# =====================================================================================
#
# Description:
#   This script performs a paired edgeR likelihood-ratio differential
#   expression analysis comparing pre- vs. posttreatment gene expression,
#   restricted to patients whose tumors were classified as PAM50 Basal
#   subtype both before and after treatment. Each patient contributes a
#   matched pre/post sample pair, and the model (~Patient + Tissue) accounts
#   for the paired (within-patient) design when testing the Pre vs. Post
#   effect.
#
# Input data:
#   mini.meta: clinical/molecular metadata data frame. Required columns:
#                 Patient_Num             - unique patient identifier
#                 PAM50_Pre_with_Clow     - PAM50 subtype, pretreatment
#                 PAM50_Post_with_Clow    - PAM50 subtype, posttreatment
#                 RNAID_Pre               - RNA-seq sample ID, pretreatment
#                 RNAID_Post              - RNA-seq sample ID, posttreatment
#
#   data.gene.raw:  raw gene expression count matrix (genes x samples), with
#                   column names matching mini.meta's RNAID_Pre/RNAID_Post
#                   sample IDs.
#
#   Gene expression data and clinical metadata are available through dbGaP
#   under accession phs003801.
#
# Output:
#   results.df: data frame of differential expression results (Post vs. Pre,
#               within Basal-Basal patients), as returned by edgeR::topTags().
#
# ==============================================================================

# ---- 1. Load required packages ----------------------------------------------
library(limma)
library(edgeR)

# ---- 2. Build the pre/post sample-pair table ----------------------------------
# Restrict to patients who were PAM50 Basal subtype both pre- and post-treatment.
basal_meta <- mini.meta[
  mini.meta$PAM50_Pre_with_Clow %in% "Basal" & mini.meta$PAM50_Post_with_Clow %in% "Basal",
  c("Patient_Num", "RNAID_Pre", "RNAID_Post")
]
n_patients <- nrow(basal_meta)

# Each patient contributes two rows: their pretreatment sample and their
# posttreatment sample, in that order.
basal.pairs <- data.frame(
  Sample = c(basal_meta$RNAID_Pre, basal_meta$RNAID_Post),
  Patient_ID = rep(basal_meta$Patient_Num, times = 2),
  PrePost = rep(c("Pre", "Post"), each = n_patients),
  stringsAsFactors = FALSE
)

# Paired index: each patient's Pre and Post sample share the same Index value
# (required by the ~Patient + Tissue paired design below).
basal.pairs$Index <- match(basal.pairs$Patient_ID, unique(basal.pairs$Patient_ID))

# ---- 3. Subset and filter the count matrix ------------------------------------
data.gene.raw.basal.pp <- data.gene.raw[, colnames(data.gene.raw) %in% basal.pairs$Sample]

# Drop versioned/ambiguous gene IDs (e.g., "ENSG00000012048.20"-style suffixes).
data.gene.raw.basal.pp <- data.gene.raw.basal.pp[!grepl("\\.[0-9]+", rownames(data.gene.raw.basal.pp)), ]

# ---- 4. Build the edgeR object and filter low-expression genes ---------------
y <- DGEList(counts = data.gene.raw.basal.pp, remove.zeros = TRUE)

keep <- suppressWarnings(filterByExpr(y))
y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)

# ---- 5. Align sample metadata to the (possibly reordered/filtered) matrix ----
basal.pairs <- basal.pairs[match(colnames(y), basal.pairs$Sample), ]

Patient <- factor(basal.pairs$Index)
Tissue <- factor(basal.pairs$PrePost, levels = c("Pre", "Post"))

# ---- 6. Estimate dispersion and fit the paired model --------------------------
design <- model.matrix(~ Patient + Tissue)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design = design, robust = TRUE)

fit <- glmFit(y, design = design)
lrt <- glmLRT(fit)  # tests the last model coefficient: TissuePost vs. TissuePre

# ---- 7. Extract results --------------------------------------------------------
results.df <- topTags(lrt, n = Inf)$table
