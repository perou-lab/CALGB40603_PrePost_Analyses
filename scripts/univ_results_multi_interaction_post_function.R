#!/usr/bin/env Rscript
# ==============================================================================
# univ_results_multi_interaction_post_function.R
# Helper function for extracting Cox model results — module x RCBCLASS interaction
# ==============================================================================
#
# Description:
#   Defines univ_results_multi_interaction_post(), which converts a list of
#   fitted coxph() model objects (survival package) into a list of tidy,
#   per-model result vectors. It is the interaction-model counterpart to
#   univ_results_multi() (see univ_results_function.R): where that function
#   summarizes models with an additive adjustment covariate
#   (~ module_score + RCBCLASS), this one summarizes models where the
#   module score's effect is allowed to differ by RCBCLASS
#   (~ module_score * RCBCLASS). See Figure5_SurvivalAnalyses.R for how
#   this is used (Analysis 4: RCB as an interaction term).
#
#   This function expects the covariate of interest (e.g., module score) to
#   be interacted with a factor that has EXACTLY 3 levels (as in the
#   70-sample post-treatment RD cohort, where RCBCLASS has levels
#   RCB-I/RCB-II/RCB-III, RCB-I as reference), producing one main-effect
#   coefficient for the covariate and two module:RCBCLASS interaction
#   coefficients. If the expected terms aren't found, it
#   stops with a clear message rather than silently returning results from
#   the wrong row.
#
#   IMPORTANT interpretation note: the "RCB2"/"RCB3" beta, HR, and
#   confidence interval fields below are the module:RCBCLASS INTERACTION
#   coefficients — i.e., the DIFFERENCE in the module score's effect in
#   that RCBCLASS group relative to the reference group (RCB-I) — not the
#   module score's total/absolute effect within that group.
#
# Parameters:
#   models : a list of fitted coxph() model objects, each fit with the
#            covariate of interest interacted with a 3-level RCBCLASS
#            factor, and with that covariate listed FIRST in the formula
#            (e.g., ~ module_score * RCBCLASS). Can be a named list; names
#            are carried through to the returned list.
#
# Returns:
#   A list (one element per model) of named numeric vectors:
#     beta_RCB1, HR_RCB1, HR.lower.confint_RCB1, HR.upper.confint_RCB1, p.value_RCB1,
#     beta_RCB2, HR_RCB2, HR.lower.confint_RCB2, HR.upper.confint_RCB2, p.value_RCB2,
#     beta_RCB3, HR_RCB3, HR.lower.confint_RCB3, HR.upper.confint_RCB3, p.value_RCB3,
#     wald.test, wald.p.value, logrank.p.value
#   (wald.test / wald.p.value / logrank.p.value describe the overall model;
#   see the interpretation note above for what the RCB1/RCB2/RCB3 fields mean.)
#
# Dependencies:
#   Requires the 'survival' package to be installed, since this function
#   calls summary() on coxph objects (an S3 method defined by 'survival').
#   Fitting the models themselves (elsewhere, before calling this function)
#   already requires 'survival' to be attached, so it isn't re-attached
#   here — but its availability is checked below so this script fails with
#   a clear message rather than a cryptic error if it's missing entirely.
#
# Sourcing:
#   This script only defines a function — it does not read or write any
#   files or reference a working directory — so it can be source()'d from
#   any directory.
#
# ==============================================================================

if (!requireNamespace("survival", quietly = TRUE)) {
  stop(
    "The 'survival' package is required by univ_results_multi_interaction_post() ",
    "(it calls summary() on fitted coxph() model objects), but it is not installed.\n",
    "Install it with: install.packages(\"survival\")"
  )
}

# ---- univ_results_multi_interaction_post(): summarize module x RCBCLASS models ---
#' Summarize a list of coxph() models with a module x RCBCLASS interaction
#'
#' @param models A list of fitted coxph() models, each fit as
#'   ~ module_score * RCBCLASS (covariate of interest listed first) where
#'   RCBCLASS is a factor with exactly 3 levels. See the file-level
#'   documentation above for what the RCB1/RCB2/RCB3 fields represent.
#' @return A list (one element per model) of named numeric vectors — see
#'   the file-level documentation above for the full field list.
univ_results_multi_interaction_post <- function(models = NA) {
  lapply(models, function(model) {
    x <- summary(model)
    coef_names <- rownames(x$coefficients)
    
    # The covariate of interest is the first term in the model's formula
    # (this function's assumed structure: covariate * RCBCLASS).
    covariate <- attr(model$terms, "term.labels")[1]
    
    main_row <- which(coef_names == covariate)
    interaction_rows <- grep(paste0("^", covariate, ":"), coef_names)
    
    if (length(main_row) != 1) {
      stop(
        "univ_results_multi_interaction_post(): expected exactly one main-",
        "effect coefficient row named '", covariate, "', found ",
        length(main_row), ". Check that the covariate of interest is the ",
        "first term in the model formula."
      )
    }
    if (length(interaction_rows) != 2) {
      stop(
        "univ_results_multi_interaction_post() expects RCBCLASS to have ",
        "exactly 3 levels (i.e., 2 '", covariate, ":RCBCLASS...' interaction ",
        "terms). Found ", length(interaction_rows), " for this model -- if ",
        "RCBCLASS has a different number of levels (e.g., a cohort that ",
        "includes \"RCB-0\"), this function's RCB1/RCB2/RCB3 output layout ",
        "does not apply as-is."
      )
    }
    
    # RCB1 = covariate's main effect (== its effect within the reference
    # RCB-I group); RCB2/RCB3 = the two interaction terms, in the order
    # RCBCLASS's factor levels were defined.
    rows <- c(main_row, interaction_rows)
    
    beta <- signif(x$coefficients[rows, 1], digits = 4)
    HR <- signif(x$coefficients[rows, 2], digits = 4)
    HR.lower <- signif(x$conf.int[rows, "lower .95"], 4)
    HR.upper <- signif(x$conf.int[rows, "upper .95"], 4)
    p.value <- signif(x$coefficients[rows, 5], digits = 4)
    
    wald.test <- signif(x$waldtest["test"], digits = 4)
    wald.pvalue <- signif(x$waldtest["pvalue"], digits = 4)
    log.p.value <- signif(x$logtest["pvalue"], digits = 4)
    
    stat_names <- c("beta", "HR", "HR.lower.confint", "HR.upper.confint", "p.value")
    group_names <- c("RCB1", "RCB2", "RCB3")
    
    # rbind() + c() interleaves per group (all 5 stats for RCB1, then RCB2,
    # then RCB3)
    res <- c(rbind(beta, HR, HR.lower, HR.upper, p.value), wald.test, wald.pvalue, log.p.value)
    names(res) <- c(
      paste0(rep(stat_names, times = length(group_names)), "_", rep(group_names, each = length(stat_names))),
      "wald.test", "wald.p.value", "logrank.p.value"
    )
    res
  })
}
