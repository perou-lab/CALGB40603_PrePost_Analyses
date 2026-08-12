#!/usr/bin/env Rscript
# ==============================================================================
# Figure 5: Cox Proportional Hazards Survival Analyses — Gene Signature Modules
# ==============================================================================
#
# Description:
#   This script runs four univariate/multivariable Cox survival analyses of
#   gene signature (module) scores against event-free survival (EFS) and
#   overall survival (OS):
#     1. All pretreatment samples (N = 340), unadjusted
#     2. Posttreatment residual disease samples (N = 70), unadjusted
#     3. Posttreatment residual disease samples (N = 70), RCB-adjusted
#     4. Posttreatment residual disease samples (N = 70), RCB as an
#        interaction term
#   Each module is tested independently (one Cox model per module), with the
#   shared covariate(s) noted above added to the adjusted/interaction models.
#
# Input data:
#   survival.meta: survival/clinical metadata data frame. Required columns:
#                     Sample                    - sample ID
#                     PreAll                    - "Yes"/"No", membership in the all-pre-treatment (N=340) cohort
#                     Post70                    - "Yes"/"No", membership in the post-treatment RD (N=70) cohort
#                     efs_time_days, efs_stat   - event-free survival time/status
#                     os_time_days, os_stat     - overall survival time/status
#
#   mini.meta: clinical metadata data frame. Required columns:
#                     RNAID_Pre, RNAID_Post     - RNA-seq sample IDs
#                     RCBCLASS                  - Residual Cancer Burden class
#
#   results.modules.published:  matrix/data frame of gene signature (module)
#                               scores, with modules as rows and sample IDs
#                               (matching survival.meta$Sample) as column names.
#
#   Gene expression data and clinical metadata are available through dbGaP
#   under accession phs003801.
#
# Sourced helper functions (assumed to live alongside this script):
#   univ_results_function.R                    
#         defines univ_results() and univ_results_multi(), which extract 
#         per-module Cox model statistics (hazard ratio, Wald p-value,
#         log-rank p-value) from a list of fitted coxph() models into a results
#         table. See that script for details.

#   univ_results_multi_interaction_post_function.R 
#         defines univ_results_multi_interaction_post(), the equivalent extraction
#         function for models with a module x RCBCLASS interaction term. See
#         that script for details.
#
# Output:
#   Eight tab-delimited Cox results tables (EFS + OS, for all four analyses),
#   written under output_dir (set below) — update this to your own output
#   location.
#
# ==============================================================================

# ---- 1. Load required packages ----------------------------------------------
library(survival)

# ---- 2. Source helper functions (see descriptions above) ---------------------
source("univ_results_function.R")
source("univ_results_multi_interaction_post_function.R")

# ---- 3. Output location -------------------------------------------------------
# Generic placeholder — set this to wherever results should be written.
output_dir <- "./Output/"

# ---- 4. Local helper: prepare a module-score matrix for Cox modeling ---------
# Subsets results.modules.published to the given samples, transposes it to
# one row per sample, and sanitizes module names so they're valid as model
# formula terms (e.g., a module named "1. Immune=Response" -> "Module_1._Immune_Response").
prepare_module_scores <- function(module_matrix, sample_ids) {
  x <- as.data.frame(t(module_matrix[, colnames(module_matrix) %in% sample_ids]))
  colnames(x) <- gsub(" ", "_", colnames(x))
  colnames(x) <- gsub("=", "", colnames(x))
  colnames(x) <- gsub("^(\\d)", "Module_\\1", colnames(x))
  colnames(x) <- gsub("/", "_", colnames(x))
  colnames(x) <- gsub("-", "_", colnames(x))
  colnames(x) <- gsub("$", "_", colnames(x), fixed = TRUE)
  x
}

# ==============================================================================
# Data preparation
# ==============================================================================
# Modifications shared across multiple analyses below are done once here,
# rather than being repeated in each analysis block.

# ---- Cohort definitions --------------------------------------------------------
survival.meta.allpre <- survival.meta[survival.meta$PreAll == "Yes", ]
survival.meta.70.post <- survival.meta[survival.meta$Post70 == "Yes", ]

# RCB-adjusted variant of the post-treatment RD cohort (used by Analyses 3 and 4).
# Computed once here since both analyses used identical input.
survival.meta.70.post.rcb <- survival.meta.70.post
survival.meta.70.post.rcb$RCBCLASS <- factor(
  mini.meta$RCBCLASS[match(survival.meta.70.post.rcb$Sample, mini.meta$RNAID_Post)],
  levels = c("RCB-I", "RCB-II", "RCB-III")
)

# ---- Module score matrices ------------------------------------------------------
module_scores_allpre <- prepare_module_scores(results.modules.published, survival.meta.allpre$Sample)
module_scores_70post <- prepare_module_scores(results.modules.published, survival.meta.70.post$Sample)

# The set of modules (and therefore covariate names) is the same across both
# cohorts — it's determined by the rows of results.modules.published, not by
# which samples/columns were selected — so a single covariates vector is used
# for every analysis below.
covariates <- colnames(module_scores_allpre)


# ---- Analysis-ready data frames (module scores + matched survival metadata) --
data_allpre <- cbind(
  module_scores_allpre,
  survival.meta.allpre[match(rownames(module_scores_allpre), survival.meta.allpre$Sample), ]
)
data_70post <- cbind(
  module_scores_70post,
  survival.meta.70.post[match(rownames(module_scores_70post), survival.meta.70.post$Sample), ]
)
data_70post_rcb <- cbind(
  module_scores_70post,
  survival.meta.70.post.rcb[match(rownames(module_scores_70post), survival.meta.70.post.rcb$Sample), ]
)

# ==============================================================================
# Analysis 1: All Pre-Treatment (N = 340), Unadjusted
# ==============================================================================

univ_formulas_efs <- sapply(covariates, function(x) as.formula(paste("Surv(efs_time_days, efs_stat)~", x)))
univ_formulas_os <- sapply(covariates, function(x) as.formula(paste("Surv(os_time_days, os_stat)~", x)))

univ_models_efs <- lapply(univ_formulas_efs, function(x) coxph(x, data = data_allpre))
univ_models_os <- lapply(univ_formulas_os, function(x) coxph(x, data = data_allpre))

univ_results_efs_preall_unadj <- as.data.frame(t(as.data.frame(univ_results(univ_models_efs), check.names = FALSE)))
univ_results_efs_preall_unadj$FDRwald <- p.adjust(univ_results_efs_preall_unadj$p.value, "fdr")
univ_results_efs_preall_unadj$FDRlogrank <- p.adjust(univ_results_efs_preall_unadj$logrank.p.value, "fdr")

univ_results_os_preall_unadj <- as.data.frame(t(as.data.frame(univ_results(univ_models_os), check.names = FALSE)))
univ_results_os_preall_unadj$FDRwald <- p.adjust(univ_results_os_preall_unadj$p.value, "fdr")
univ_results_os_preall_unadj$FDRlogrank <- p.adjust(univ_results_os_preall_unadj$logrank.p.value, "fdr")

write.table(
  univ_results_efs_preall_unadj,
  file.path(output_dir, "EFS_Coxph_Results_All_Pre_Unadjusted.txt"),
  sep = "\t", col.names = NA
)

write.table(
  univ_results_os_preall_unadj,
  file.path(output_dir, "OS_Coxph_Results_All_Pre_Unadjusted.txt"),
  sep = "\t", col.names = NA
)

# ==============================================================================
# Analysis 2: Post-Treatment Residual Disease (N = 70), Unadjusted
# ==============================================================================

univ_formulas_efs <- sapply(covariates, function(x) as.formula(paste("Surv(efs_time_days, efs_stat)~", x)))
univ_formulas_os <- sapply(covariates, function(x) as.formula(paste("Surv(os_time_days, os_stat)~", x)))

univ_models_efs <- lapply(univ_formulas_efs, function(x) coxph(x, data = data_70post))
univ_models_os <- lapply(univ_formulas_os, function(x) coxph(x, data = data_70post))

univ_results_efs_70post_unadj <- as.data.frame(t(as.data.frame(univ_results(univ_models_efs), check.names = FALSE)))
univ_results_efs_70post_unadj$FDRwald <- p.adjust(univ_results_efs_70post_unadj$p.value, "fdr")
univ_results_efs_70post_unadj$FDRlogrank <- p.adjust(univ_results_efs_70post_unadj$logrank.p.value, "fdr")

univ_results_os_70post_unadj <- as.data.frame(t(as.data.frame(univ_results(univ_models_os), check.names = FALSE)))
univ_results_os_70post_unadj$FDRwald <- p.adjust(univ_results_os_70post_unadj$p.value, "fdr")
univ_results_os_70post_unadj$FDRlogrank <- p.adjust(univ_results_os_70post_unadj$logrank.p.value, "fdr")

write.table(
  univ_results_efs_70post_unadj,
  file.path(output_dir, "EFS_Coxph_Results_Post_70RD_Unadjusted.txt"),
  sep = "\t", col.names = NA
)

write.table(
  univ_results_os_70post_unadj,
  file.path(output_dir, "OS_Coxph_Results_Post_70RD_Unadjusted.txt"),
  sep = "\t", col.names = NA
)

# ==============================================================================
# Analysis 3: Post-Treatment Residual Disease (N = 70), RCB-Adjusted
# ==============================================================================

univ_formulas_efs <- sapply(covariates, function(x) as.formula(paste("Surv(efs_time_days, efs_stat)~", x, "+ RCBCLASS")))
univ_formulas_os <- sapply(covariates, function(x) as.formula(paste("Surv(os_time_days, os_stat)~", x, "+ RCBCLASS")))

univ_models_efs <- lapply(univ_formulas_efs, function(x) coxph(x, data = data_70post_rcb))
univ_models_os <- lapply(univ_formulas_os, function(x) coxph(x, data = data_70post_rcb))

univ_results_efs_70post_rcb <- as.data.frame(t(as.data.frame(univ_results_multi(univ_models_efs), check.names = FALSE)))
univ_results_efs_70post_rcb$FDR.adj.pval <- p.adjust(univ_results_efs_70post_rcb$variable.adjusted.p.value, "fdr")
univ_results_efs_70post_rcb$FDRlogrank <- p.adjust(univ_results_efs_70post_rcb$logrank.p.value, "fdr")

univ_results_os_70post_rcb <- as.data.frame(t(as.data.frame(univ_results_multi(univ_models_os), check.names = FALSE)))
univ_results_os_70post_rcb$FDR.adj.pval <- p.adjust(univ_results_os_70post_rcb$variable.adjusted.p.value, "fdr")
univ_results_os_70post_rcb$FDRlogrank <- p.adjust(univ_results_os_70post_rcb$logrank.p.value, "fdr")

write.table(
  univ_results_efs_70post_rcb,
  file.path(output_dir, "EFS_Coxph_Results_Post_70RD_RCB_adjusted.txt"),
  sep = "\t", col.names = NA
)

write.table(
  univ_results_os_70post_rcb,
  file.path(output_dir, "OS_Coxph_Results_Post_70RD_RCB_adjusted.txt"),
  sep = "\t", col.names = NA
)

# ==============================================================================
# Analysis 4: Post-Treatment Residual Disease (N = 70), RCB as Interaction Term
# ==============================================================================

univ_formulas_efs <- sapply(covariates, function(x) as.formula(paste("Surv(efs_time_days, efs_stat)~", x, "* RCBCLASS")))
univ_formulas_os <- sapply(covariates, function(x) as.formula(paste("Surv(os_time_days, os_stat)~", x, "* RCBCLASS")))

univ_models_efs <- lapply(univ_formulas_efs, function(x) coxph(x, data = data_70post_rcb))
univ_models_os <- lapply(univ_formulas_os, function(x) coxph(x, data = data_70post_rcb))

univ_results_efs_70post_interaction <- as.data.frame(t(as.data.frame(univ_results_multi_interaction_post(univ_models_efs), check.names = FALSE)))
univ_results_efs_70post_interaction$FDRlogrank <- p.adjust(univ_results_efs_70post_interaction$logrank.p.value, "fdr")

univ_results_os_70post_interaction <- as.data.frame(t(as.data.frame(univ_results_multi_interaction_post(univ_models_os), check.names = FALSE)))
univ_results_os_70post_interaction$FDRlogrank <- p.adjust(univ_results_os_70post_interaction$logrank.p.value, "fdr")

write.table(
  univ_results_efs_70post_interaction,
  file.path(output_dir, "EFS_Coxph_Results_Post_70RD_RCB_as_interaction.txt"),
  sep = "\t", col.names = NA
)

write.table(
  univ_results_os_70post_interaction,
  file.path(output_dir, "OS_Coxph_Results_Post_70RD_RCB_as_interaction.txt"),
  sep = "\t", col.names = NA
)


