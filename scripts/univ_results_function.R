#!/usr/bin/env Rscript
# ==============================================================================
# univ_results_function.R
# Helper functions for extracting Cox proportional hazards model results
# ==============================================================================
#
# Description:
#   Defines two helper functions that convert a list of fitted coxph()
#   model objects (survival package) into a list of tidy, per-model result
#   vectors — hazard ratio, confidence interval, and the relevant
#   significance tests. Intended for loop-based per-module Cox survival
#   analyses (see Figure5_SurvivalAnalyses.R), where one coxph() model is
#   fit per gene signature/module and these functions summarize all of
#   them in one pass.
#
#   univ_results(models)
#     For UNIVARIATE models — a single covariate of interest and no
#     adjustment covariates, e.g. Surv(time, status) ~ module_score.
#
#   univ_results_multi(models)
#     For MULTIVARIABLE models where the covariate of interest is listed
#     FIRST in the model formula, followed by one or more adjustment
#     covariates, e.g. Surv(time, status) ~ module_score + RCBCLASS.
#     Reports results for the covariate of interest (the first row of the
#     model's coefficient table) alongside the model's overall Wald and
#     log-rank test statistics.
#
# Parameters (both functions):
#   models : a list of fitted coxph() model objects (e.g., the output of
#            lapply(formulas, function(f) coxph(f, data = ...))). Can be a
#            named list; names are carried through to the returned list.
#
# Returns (both functions):
#   A list (one element per input model) of named numeric vectors — see
#   each function's documentation below for its exact fields.
#
# Dependencies:
#   Requires the 'survival' package to be installed, since these functions
#   call summary() on coxph objects (an S3 method defined by 'survival').
#   Fitting the models themselves (elsewhere, before calling these
#   functions) already requires 'survival' to be attached, so it isn't
#   re-attached here — but its availability is checked below so this
#   script fails with a clear message rather than a cryptic error if it's
#   missing entirely.
#
# Sourcing:
#   This script only defines functions — it does not read or write any
#   files or reference a working directory — so it can be source()'d from
#   any directory.
#
# ==============================================================================

if (!requireNamespace("survival", quietly = TRUE)) {
  stop(
    "The 'survival' package is required by univ_results() / univ_results_multi() ",
    "(they call summary() on fitted coxph() model objects), but it is not installed.\n",
    "Install it with: install.packages(\"survival\")"
  )
}

# ---- univ_results(): summarize univariate Cox models --------------------------
#' Summarize a list of univariate coxph() models
#'
#' @param models A list of fitted coxph() models, each with exactly ONE
#'   covariate (no adjustment covariates). For a multivariable model, use
#'   univ_results_multi() instead — this function assumes a single-row
#'   coefficient table and will stop with an error if that assumption
#'   doesn't hold, rather than silently returning results for the wrong
#'   covariate.
#' @return A list (one element per model) of named numeric vectors:
#'   beta, HR, HR.lower.confint, HR.upper.confint, wald.test, p.value,
#'   logrank.p.value.
univ_results <- function(models = NA) {
  lapply(models, function(x) {
    x <- summary(x)
    
    if (nrow(x$coefficients) != 1) {
      stop(
        "univ_results() expects models with exactly one covariate (a ",
        "univariate model). This model has ", nrow(x$coefficients),
        " covariates -- use univ_results_multi() for adjusted/",
        "multivariable models instead."
      )
    }
    
    p.value <- signif(x$waldtest["pvalue"], digits = 3)
    log.p.value <- signif(x$logtest["pvalue"], digits = 3)
    wald.test <- signif(x$waldtest["test"], digits = 3)
    beta <- signif(x$coefficients[1], digits = 3)
    HR <- signif(x$coefficients[2], digits = 4)
    HR.confint.lower <- signif(x$conf.int[, "lower .95"], 4)
    HR.confint.upper <- signif(x$conf.int[, "upper .95"], 4)
    
    res <- c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, p.value, log.p.value)
    names(res) <- c(
      "beta", "HR", "HR.lower.confint", "HR.upper.confint",
      "wald.test", "p.value", "logrank.p.value"
    )
    res
  })
}

# ---- univ_results_multi(): summarize multivariable Cox models -----------------
#' Summarize a list of multivariable coxph() models
#'
#' @param models A list of fitted coxph() models, each with the covariate
#'   of interest as the FIRST term in the model formula, followed by one
#'   or more adjustment covariates (e.g., ~ module_score + RCBCLASS).
#' @return A list (one element per model) of named numeric vectors:
#'   beta, HR, HR.lower.confint, HR.upper.confint, wald.test, wald.p.value,
#'   logrank.p.value, variable.adjusted.p.value. beta / HR / confidence
#'   interval / variable.adjusted.p.value describe the first (primary)
#'   covariate only; wald.test / wald.p.value / logrank.p.value describe
#'   the overall model.
univ_results_multi <- function(models = NA) {
  lapply(models, function(x) {
    x <- summary(x)
    
    p.value <- signif(x$coefficients[1, 5], digits = 3)
    log.p.value <- signif(x$logtest["pvalue"], digits = 3)
    wald.pvalue <- signif(x$waldtest["pvalue"], digits = 3)
    wald.test <- signif(x$waldtest["test"], digits = 3)
    beta <- signif(x$coefficients[1, 1], digits = 3)
    HR <- signif(x$coefficients[1, 2], digits = 4)
    HR.confint.lower <- signif(x$conf.int[1, "lower .95"], 4)
    HR.confint.upper <- signif(x$conf.int[1, "upper .95"], 4)
    
    res <- c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, wald.pvalue, log.p.value, p.value)
    names(res) <- c(
      "beta", "HR", "HR.lower.confint", "HR.upper.confint", "wald.test",
      "wald.p.value", "logrank.p.value", "variable.adjusted.p.value"
    )
    res
  })
}
