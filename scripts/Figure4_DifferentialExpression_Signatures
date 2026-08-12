#!/usr/bin/env Rscript
# ==================================================================================================
# Figure 4 (Signatures): SAM Input Preparation — Paired Pre- vs. Posttreatment,
#                         Patients Who Remained Basal Subtype Throughout: Basal-Basal (N = 21 pairs) 
# ==================================================================================================
#
# Description:
#   This script prepares a SAM (Significance Analysis of Microarrays) PAIRED
#   input file comparing pre- vs. posttreatment gene signature/module scores,
#   restricted to patients whose tumors were classified as PAM50 Basal
#   subtype both before and after treatment. Each patient's pretreatment
#   sample is labeled -i and their post-treatment sample is labeled +i (SAM's
#   paired-data response format), so that SAM analyzes each patient as a
#   matched before/after pair rather than as two independent groups.
#
#   The script writes a tab-delimited file that must be manually reformatted
#   in Excel (see instructions below) before running SAM. runSAM() is part of 
#   the samr package. 
#
# Input data:
#   mini.meta: clinical/molecular metadata data frame. Required columns:
#                 RNAID_Pre               - RNA-seq sample ID, pre-treatment
#                 RNAID_Post              - RNA-seq sample ID, post-treatment ("." if unavailable)
#                 PAM50_Pre_with_Clow     - PAM50 subtype, pre-treatment
#                 PAM50_Post_with_Clow    - PAM50 subtype, post-treatment
#
#   module.results.published: matrix/data frame of gene signature (module)
#                             scores, with signatures as rows and sample IDs
#                             (matching RNAID_Pre / RNAID_Post) as column names.
#
#   Gene expression data and clinical metadata are available through dbGaP
#   under accession phs003801.
#
# Output:
#   Modules_BasalBasalPre-Clow_vs_Post_SAMR.txt : tab-delimited signature score
#     matrix with sample columns labeled by SAM paired response code (-i for
#     patient i's pre-treatment sample, +i for their post-treatment sample),
#     written to:
#     ./Output/File/Path/
#
# ==============================================================================

# ---- 1. Identify the paired Pre/Post samples -----------------------------------
# Restrict to patients who were PAM50 Basal subtype both pre- and post-treatment.
basal_meta <- mini.meta[
  !(mini.meta$RNAID_Post %in% ".") &
    mini.meta$PAM50_Pre_with_Clow %in% "Basal" &
    mini.meta$PAM50_Post_with_Clow %in% "Basal",
  c("RNAID_Pre", "RNAID_Post")
]
n_patients <- nrow(basal_meta)

# SAM paired-data response codes: 
# pretreatment sample = -i.
# posttreatment sample = +i.
group <- setNames(
  c(-seq_len(n_patients), seq_len(n_patients)),
  c(basal_meta$RNAID_Pre, basal_meta$RNAID_Post)
)

# ---- 2. Subset signature scores to the paired samples ---------------------------
module.results.published.prepost.sam <- module.results.published[
  , colnames(module.results.published) %in% names(group)
]

# ---- 3. Relabel sample columns with their SAM paired response code -------------
colnames(module.results.published.prepost.sam) <- group[
  match(colnames(module.results.published.prepost.sam), names(group))
]

# ---- 4. Write the SAM input file ------------------------------------------------
write.table(
  module.results.published.prepost.sam,
  "./Output/File/Path/file.txt",
  sep = "\t", col.names = NA
)

# ==============================================================================
# Required Excel modifications before running SAM
# (per the official SAM "Significance Analysis of Microarrays" User Guide,
#  Stanford: http://www-stat.stanford.edu/~tibs/SAM)
# ==============================================================================
#
# The file written above has ONE identifier column (signature name, from
# write.table(col.names = NA)) followed directly by the sample columns. SAM's
# required layout has TWO identifier columns before the data begins:
#
#   Column A  = Gene/Signature Name
#   Column B  = Gene/Signature ID
#   Column C+ = signature score values, one column per sample, with Row 1
#               containing that sample's paired response code (-i or +i)
#   Row 1     = response codes, starting at COLUMN C
#   Row 2+    = one row per signature, with the score values as plain numbers
#
# Steps to take in Excel after opening the .txt file:
#   1. Open the file in Excel (it will import as tab-delimited).
#   2. Insert a new column B, and copy column A (signature names) into it, so
#      there are two identifier columns (A: name, B: ID) before the sample data.
#      This shifts all sample columns from starting at column B to column C.
#   3. Confirm that Row 1 of the sample columns (C onward) contains only the 
#      numeric paired response codes (e.g., -1, 1, -2, 2, ...) — not any other
#      header text.
#   4. Excel may interpret a leading "-" or "+" as the start of a formula; if a
#      response code cell shows an error or blank instead of the number,
#      re-enter it as text-formatted or prefix with an apostrophe so it reads
#      as -1 rather than a formula.
#   5. Confirm all signature score values are stored as numbers.
#   6. Save the file as an Excel workbook (.xls/.xlsx) — SAM's Excel Add-in
#      operates on an open Excel worksheet, not the original .txt file.
#
# When running SAM from the Excel Add-in dialog, select "Paired data" as the
# response type (NOT "Two class, unpaired") — this script's response codes
# only make sense under the paired analysis, which tests each patient's own
# pre-to-post change rather than treating all Pre/Post samples as independent.
#
# Once the Excel worksheet is prepared, proceed to:
library(samr)
runSAM()


