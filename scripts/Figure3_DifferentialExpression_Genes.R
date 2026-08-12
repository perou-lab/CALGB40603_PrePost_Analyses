#!/usr/bin/env Rscript
# ======================================================================================================
# Figure 3: Differential Gene Expression — Pretreatment Basal-Basal (N= 21) vs. Basal-NonBasal (N = 102)
# ======================================================================================================
#
# Description:
#   This script performs an edgeR quasi-likelihood differential expression
#   analysis comparing pretreatment gene expression between two groups of
#   patients whose tumors were classified as PAM50 Basal at baseline:
#     - "BasalBasal"    : remained Basal subtype after treatment
#     - "BasalNonBasal" : switched to a different PAM50 subtype after treatment
#   The comparison uses each patient's pre-treatment ("RNAID_Pre") sample, and
#   is restricted to patients with a paired post-treatment sample available.
#
# Input data:
#   mini.meta: clinical/molecular metadata data frame. Required columns:
#                 RNAID_Pre             - RNA-seq sample ID, pretreatment
#                 RNAID_Post            - RNA-seq sample ID, posttreatment ("." if unavailable)
#                 PAM50_Pre_with_Clow   - PAM50 subtype, pretreatment
#                 PAM50_Post_with_Clow  - PAM50 subtype, posttreatment
#
#   data.gene.raw : raw gene expression count matrix (genes x samples), with
#                   column names matching mini.meta's RNAID_Pre sample IDs.
#
#   Gene expression data and clinical metadata are available through dbGaP
#   under accession phs003801.
#
# Output:
#   results.df : data frame of differential expression results (BasalBasal vs.
#                BasalNonBasal), as returned by edgeR::topTags().
#
# ==============================================================================

# ---- 1. Load required packages ----------------------------------------------
library(limma)
library(edgeR)

# ---- 2. Identify the comparison groups ---------------------------------------
# Both groups are PAM50 Basal subtype at baseline (pretreatment) with a
# paired post-treatment sample; the groups differ in whether the tumor
# remained Basal or switched to a different subtype after treatment.
basal_pre <- mini.meta[mini.meta$RNAID_Post != "." & mini.meta$PAM50_Pre_with_Clow == "Basal", ]

sample_ids <- basal_pre$RNAID_Pre
group <- factor(
  ifelse(basal_pre$PAM50_Post_with_Clow == "Basal", "BasalBasal", "BasalNonBasal"),
  levels = c("BasalNonBasal", "BasalBasal")
)

# Keep only samples actually present in the count matrix.
in_matrix <- sample_ids %in% colnames(data.gene.raw)
sample_ids <- sample_ids[in_matrix]
group <- group[in_matrix]

# ---- 3. Subset and filter the count matrix ------------------------------------
data.gene.raw.basal <- data.gene.raw[, sample_ids]

# Drop versioned/ambiguous gene IDs (e.g., "ENSG00000012048.20"-style suffixes).
data.gene.raw.basal <- data.gene.raw.basal[!grepl("\\.[0-9]+", rownames(data.gene.raw.basal)), ]

# ---- 4. Build the edgeR object and filter low-expression genes ---------------
y <- DGEList(counts = data.gene.raw.basal, group = group, remove.zeros = TRUE)

keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)

# ---- 5. Estimate dispersion and fit the quasi-likelihood model ---------------
design <- model.matrix(~group)
y <- estimateDisp(y, design = design, robust = TRUE)

fit <- glmQLFit(y, design = design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)

# ---- 6. Extract results --------------------------------------------------------
results.df <- topTags(qlf, n = Inf)$table









