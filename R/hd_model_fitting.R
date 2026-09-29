# Portions of this file are adapted from the FAIBCompiler R package
# (https://github.com/bcgov/FAIBCompiler), files R/DBH_Height_MEM.R and
# R/heightEstimate_byHeightModel.R.
# Copyright 2019 Province of British Columbia. Original author: Yong Luo.
# Licensed under the Apache License, Version 2.0; see the LICENSE file.
# Modified 2026 by Natural Resources Canada: generalised to six H-D model
# forms with user-selectable fixed- or mixed-effects fitting, automatic
# starting values, BLUP calibration, broken-top handling, and a documented
# SP0 -> SP_TYPE fallback hierarchy.

# Fitting height-diameter models to PSP data.  Exported functions:
# fit_hd_model()              - fit a single-group H-D model (nls or nlme)
# fit_hd_models_by_group()   - SP0 + SP_TYPE hierarchical fitting pipeline
# ht_impute()                 - impute missing heights from fit_hd_models_by_group() output
# hd_start_values()          - automatic starting-value computation
# ht_predict()               - apply a hierarchy of coef tables to tree data
# calibrate_hd_blup()        - BLUP calibration of a mixed-effects model to a new plot

# Internal start-value helpers -----------------------------------------------
# Each function returns a named list suitable for nls(start = ...) or
# as.vector for nlme(start = ...).

.start_naslund <- function(dbh, height) {
  # Linearise: DBH / sqrt(H - 1.3) = a + b * DBH
  y   <- dbh / sqrt(pmax(height - 1.3, 0.01))
  fit <- stats::lm(y ~ dbh)
  cf  <- unname(stats::coef(fit))
  list(a = max(cf[1], 0.01), b = max(cf[2], 1e-4))
}

.start_curtis <- function(dbh, height) {
  # log(H - 1.3) = log(a) + b * log(DBH / (1 + DBH))
  y   <- log(pmax(height - 1.3, 0.01))
  x   <- log(dbh / (1 + dbh))
  fit <- stats::lm(y ~ x)
  cf  <- unname(stats::coef(fit))
  list(a = max(exp(cf[1]), 0.01), b = max(cf[2], 1e-4))
}

.start_logistic <- function(dbh, height) {
  list(a = max(height) * 1.1 - 1.3, b = 1.0, c = 0.1)
}

.start_korf <- function(dbh, height) {
  list(a = max(height) * 1.1 - 1.3, b = 1.0, c = 0.5)
}

.start_weibull <- function(dbh, height) {
  list(a = max(height) * 1.1 - 1.3, b = 0.05, c = 1.0)
}

.start_richards <- function(dbh, height) {
  list(a = max(height) * 1.1 - 1.3, b = 0.05, c = 1.5)
}


# Supported model formulae (height as response, DBH as predictor) ------------

.hd_formula <- function(model) {
  switch(model,
    naslund  = stats::as.formula("HEIGHT ~ 1.3 + (DBH^2) / (a + b * DBH)^2"),
    curtis   = stats::as.formula("HEIGHT ~ 1.3 + a * (DBH / (1 + DBH))^b"),
    logistic = stats::as.formula("HEIGHT ~ 1.3 + a / (1 + b * exp(-c * DBH))"),
    korf     = stats::as.formula("HEIGHT ~ 1.3 + a * exp(-b * DBH^(-c))"),
    weibull  = stats::as.formula("HEIGHT ~ 1.3 + a * (1 - exp(-b * DBH^c))"),
    richards = stats::as.formula("HEIGHT ~ 1.3 + a * (1 - exp(-b * DBH))^c")
  )
}


#' Compute starting values for height-diameter model fitting
#'
#' @description
#' Returns a named list of reasonable starting parameter values for the
#' specified height-diameter model form. Useful when calling \code{nls} or
#' \code{nlme} directly. Two-parameter models (\code{naslund}, \code{curtis})
#' use a linearisation strategy; three-parameter models use heuristics based
#' on the observed height range.
#'
#' @param dbh    numeric vector. Diameter at breast height (cm).
#' @param height numeric vector. Total tree height (m). Must be the same length
#'   as \code{dbh}.
#' @param model  character. One of \code{"naslund"}, \code{"curtis"},
#'   \code{"logistic"}, \code{"korf"}, \code{"weibull"}, \code{"richards"}.
#'   Case-insensitive.
#'
#' @return A named list of starting values, e.g. \code{list(a = ..., b = ...)}
#'   for two-parameter models or \code{list(a = ..., b = ..., c = ...)} for
#'   three-parameter models.
#' @export
#' @seealso \code{\link{fit_hd_model}}
#' @examples
#' set.seed(1)
#' dbh <- runif(100, 5, 60)
#' ht  <- ht_korf(dbh, a = 38, b = 2.5, c = 0.4) + rnorm(100, 0, 1)
#' hd_start_values(dbh, ht, model = "korf")
hd_start_values <- function(dbh, height, model) {
  model <- tolower(model)
  fn <- switch(model,
    naslund  = .start_naslund,
    curtis   = .start_curtis,
    logistic = .start_logistic,
    korf     = .start_korf,
    weibull  = .start_weibull,
    richards = .start_richards,
    stop(
      "Unknown model '", model, "'. ",
      "Must be one of: naslund, curtis, logistic, korf, weibull, richards."
    )
  )
  fn(dbh, height)
}


#' Fit a height-diameter allometric model
#'
#' @description
#' Fits a nonlinear height-diameter model to a \strong{single} group of trees.
#' This is the low-level workhorse function. For the standard PSP workflow
#' (SP0-level fitting with SP_TYPE fallback and automatic model selection) use
#' \code{\link{fit_hd_models_by_group}}, which calls this function internally
#' and provides greater convenience. Use \code{fit_hd_model} directly when you
#' need precise control over a single fit (custom grouping, non-standard
#' thresholds, or step-by-step exploration).
#'
#' Starting values are computed automatically via
#' \code{\link{hd_start_values}}.
#'
#' \strong{Fixed-effects mode} (default): uses \code{\link[stats]{nls}} from
#' base R. No additional packages required.
#'
#' \strong{Mixed-effects mode}: activated by supplying \code{group_col}. A
#' site-level random effect is placed on the asymptote parameter \eqn{a} using
#' \code{nlme::nlme}. Requires the \pkg{nlme} package to be installed.
#'
#' All six supported model forms use a breast-height stump convention of 1.3 m.
#' Heights are expected in metres and DBH in centimetres.
#'
#' @param data        data frame (or data.table). Must contain columns for DBH
#'   and height (see \code{dbh_col} and \code{height_col}). Rows with
#'   \code{NA} in any required column, DBH \eqn{\leq 0}, or
#'   height \eqn{\leq 1.3} m are silently dropped before fitting.
#'   \strong{Broken-top trees must be excluded before calling this function.}
#'   A broken-top tree's recorded height is the broken-top height, not total
#'   height; including it biases the model. Filter with
#'   \code{data[!data$BTOP, ]} (or equivalent) before passing to
#'   \code{fit_hd_model}. Broken-top trees have no valid total height and
#'   must be excluded before fitting. Filter with
#'   \code{data[!data$BTOP, ]} (or equivalent) before passing to
#'   \code{fit_hd_model}.
#' @param model       character. Model form to fit. One of \code{"naslund"},
#'   \code{"curtis"}, \code{"logistic"}, \code{"korf"}, \code{"weibull"},
#'   \code{"richards"}. Default \code{"naslund"}.
#' @param height_col  character. Name of the height column. Default
#'   \code{"HEIGHT"}.
#' @param dbh_col     character. Name of the DBH column. Default \code{"DBH"}.
#' @param group_col   character or \code{NULL}. Name of a grouping column
#'   (e.g. plot, site, or species). When supplied, a mixed-effects model is
#'   fitted placing a random effect on \eqn{a} for each level of this group
#'   (requires \pkg{nlme}). Can be any column that identifies a group --
#'   plot ID, BEC zone, species, stand, etc. Default \code{NULL}.
#' @param nested_col  character or \code{NULL}. Name of a second grouping
#'   column nested within \code{group_col} (e.g. tree ID for PSP repeated
#'   measures). When supplied together with \code{group_col}, a two-level
#'   nested random effect is placed on \eqn{a}
#'   (\code{a ~ 1 | group / nested}, i.e. outer/inner in \pkg{nlme}
#'   notation). Useful for PSP data where the same tree is measured across
#'   multiple measurement periods. Unbalanced data are handled gracefully:
#'   trees measured only once receive a tree-level BLUP of zero (the model
#'   falls back to the plot-level prediction for them), so you do not need to
#'   omit \code{nested_col} just because some trees have a single measurement.
#'   \strong{True nesting is required}: the grouping structure must be
#'   hierarchical, not crossed. Every unit in \code{nested_col} must belong
#'   to exactly one level of \code{group_col}. \pkg{nlme} automatically
#'   handles non-globally-unique inner IDs (e.g. tree \dQuote{1} in plot
#'   \dQuote{A} and tree \dQuote{1} in plot \dQuote{B} are kept distinct
#'   via compound keys). If your design is genuinely crossed rather than
#'   nested (e.g. the same observer works across multiple plots), a
#'   different random-effects specification is required.
#'   Requires \code{group_col}. Default \code{NULL}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{model}}{character. The model form that was fitted.}
#'     \item{\code{fit}}{The fitted model object (\code{nls} or \code{nlme}).}
#'     \item{\code{coefficients}}{Named numeric vector of fixed-effect
#'       coefficients.}
#'     \item{\code{r2_marginal}}{Numeric. Marginal \eqn{R^2}: proportion of
#'       total height variance explained by the fixed-effects prediction.}
#'     \item{\code{r2_conditional}}{Numeric. Conditional \eqn{R^2}: proportion of
#'       total height variance explained by the full model including random
#'       effects. Computed from the innermost-level residuals (\code{TREE_ID}
#'       column when \code{nested_col} is supplied; \code{SITE_ID} column
#'       otherwise). Returned only for mixed-effects fits.}
#'     \item{\code{random_effects_outer}}{Data frame of outer-group BLUPs for
#'       \eqn{a}, one row per level of \code{group_col} (returned when
#'       \code{group_col} is supplied). \dQuote{Outer} follows \pkg{nlme}
#'       notation: the coarser, wrapping grouping factor.}
#'     \item{\code{random_effects_inner}}{Data frame of inner-group BLUPs for
#'       \eqn{a}, one row per level of \code{nested_col} (returned only when
#'       \code{nested_col} is supplied). \dQuote{Inner} is the finer grouping
#'       nested inside the outer factor.}
#'   }
#'
#' @export
#' @seealso \code{\link{ht_from_dbh}}, \code{\link{hd_start_values}},
#'   \code{\link{ht_predict}}
#' @examples
#' # Fixed-effects fit
#' set.seed(42)
#' dbh    <- runif(200, 5, 60)
#' height <- ht_naslund(dbh, a = 1.2, b = 0.04) + rnorm(200, 0, 1.5)
#' dat    <- data.frame(DBH = dbh, HEIGHT = height)
#'
#' result <- fit_hd_model(dat, model = "naslund")
#' result$coefficients
#' result$r2_marginal
#'
#' # Single-level mixed-effects fit (requires nlme)
#' \dontrun{
#' dat$PLOT <- sample(paste0("P", 1:10), 200, replace = TRUE)
#' result_me <- fit_hd_model(dat, model = "naslund", group_col = "PLOT")
#' result_me$coefficients
#' result_me$random_effects_outer
#'
#' # Two-level mixed-effects fit: group + nested (for PSP repeated measures)
#' dat$TREE_ID <- sample(paste0("T", 1:50), 200, replace = TRUE)
#' result_2l  <- fit_hd_model(dat, model = "naslund",
#'                            group_col = "PLOT", nested_col = "TREE_ID")
#' result_2l$random_effects_inner
#'
#' # -------------------------------------------------------------------
#' # Manual SP0-level loop with SP_TYPE fallback and BLUP-adjusted
#' # prediction — the pattern used by fit_hd_models_by_group() internally,
#' # exposed here for users who need full control.
#' # -------------------------------------------------------------------
#' library(CanWestAllometryCFS)
#' trees <- psp_trees
#' trees$SPECIES_CORR    <- species_correction(trees$SPECIES, trees$BEC_ZONE)
#' trees$SPECIES_SP0     <- bc_species_to_sp0(trees$SPECIES_CORR)
#' trees$SPECIES_SP_TYPE <- bc_species_to_sp_type(trees$SPECIES_CORR)
#'
#' measured <- trees[!is.na(trees$HEIGHT) & !trees$BTOP, ]
#'
#' # Helper: try ME then FE for one subset of trees
#' .fit_one <- function(sub) {
#'   n_plots <- length(unique(sub$SITE_IDENTIFIER))
#'   fit <- NULL; method <- "none"
#'   if (nrow(sub) >= 10 && n_plots >= 5 &&
#'       requireNamespace("nlme", quietly = TRUE)) {
#'     has_rep <- anyDuplicated(sub$unitreeid) > 0
#'     fit <- tryCatch(
#'       fit_hd_model(sub, model = "naslund",
#'                    group_col  = "SITE_IDENTIFIER",
#'                    nested_col = if (has_rep) "unitreeid" else NULL),
#'       error = function(e) NULL)
#'     if (!is.null(fit)) method <- "mixed-effects"
#'   }
#'   if (is.null(fit) && nrow(sub) >= 10) {
#'     fit <- tryCatch(fit_hd_model(sub, model = "naslund"),
#'                     error = function(e) NULL)
#'     if (!is.null(fit)) method <- "fixed-effects"
#'   }
#'   list(fit = fit, method = method, n = nrow(sub))
#' }
#'
#' # Level 1: one model per SP0 group
#' sp0_list <- sort(unique(measured$SPECIES_SP0))
#' hd_sp0   <- setNames(lapply(sp0_list, function(sp0)
#'   c(.fit_one(measured[measured$SPECIES_SP0 == sp0, ]),
#'     list(level = "sp0", group = sp0))),
#'   sp0_list)
#'
#' # Level 2: SP_TYPE fallback for SP0 groups with no model
#' no_model    <- vapply(hd_sp0, function(m) m$method == "none", logical(1))
#' spt_needed  <- sort(unique(
#'   measured$SPECIES_SP_TYPE[measured$SPECIES_SP0 %in% names(no_model)[no_model]]))
#' hd_sptype   <- setNames(lapply(spt_needed, function(spt)
#'   c(.fit_one(measured[measured$SPECIES_SP_TYPE == spt, ]),
#'     list(level = "sp_type", group = spt))),
#'   spt_needed)
#'
#' # Imputation using fixed + BLUP-adjusted 'a'
#' .blup_a <- function(ids, re) {
#'   if (is.null(re)) return(rep(0, length(ids)))
#'   b <- re$a[match(ids, re$SITE_IDENTIFIER)]
#'   ifelse(is.na(b), 0, b)
#' }
#' trees$HT_PROJ <- trees$HEIGHT
#' for (m in c(hd_sp0, hd_sptype)) {
#'   if (is.null(m$fit)) next
#'   cf   <- m$fit$coefficients
#'   re   <- m$fit$random_effects_site
#'   miss <- if (m$level == "sp0")
#'             is.na(trees$HEIGHT) & !trees$BTOP &
#'               !is.na(trees$SPECIES_SP0) & trees$SPECIES_SP0 == m$group
#'           else
#'             is.na(trees$HT_PROJ) & !trees$BTOP &
#'               !is.na(trees$SPECIES_SP_TYPE) & trees$SPECIES_SP_TYPE == m$group
#'   if (!any(miss)) next
#'   a_pred <- cf["a"] + .blup_a(trees$SITE_IDENTIFIER[miss], re)
#'   trees$HT_PROJ[miss] <- round(
#'     ht_from_dbh(trees$DBH[miss], "naslund", a = a_pred, b = cf["b"]), 1)
#' }
#' }
fit_hd_model <- function(data,
                         model       = "naslund",
                         height_col  = "HEIGHT",
                         dbh_col     = "DBH",
                         group_col   = NULL,
                         nested_col  = NULL) {
  model <- tolower(model)

  if (!is.null(nested_col) && is.null(group_col)) {
    stop("'nested_col' requires 'group_col' to also be supplied.")
  }

  # Validate columns exist
  required <- c(height_col, dbh_col, group_col, nested_col)
  missing  <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Column(s) not found in data: ", paste(missing, collapse = ", "))
  }

  # Filter valid rows
  keep <- !is.na(data[[dbh_col]])    &
          !is.na(data[[height_col]]) &
          data[[dbh_col]]    > 0     &
          data[[height_col]] > 1.3
  if (!is.null(group_col)) {
    keep <- keep & !is.na(data[[group_col]])
  }
  if (!is.null(nested_col)) {
    keep <- keep & !is.na(data[[nested_col]])
  }

  n_valid <- sum(keep)
  if (n_valid < 10) {
    stop("Fewer than 10 valid observations after filtering. Cannot fit model.")
  }

  # Build a clean data frame with standardised column names
  df <- data.frame(
    DBH    = data[[dbh_col]][keep],
    HEIGHT = data[[height_col]][keep]
  )
  if (!is.null(group_col)) {
    df$SITE_IDENTIFIER <- data[[group_col]][keep]
  }
  if (!is.null(nested_col)) {
    df$unitreeid <- data[[nested_col]][keep]
  }

  start        <- hd_start_values(df$DBH, df$HEIGHT, model)
  formula_fit  <- .hd_formula(model)
  R2_bottom    <- sum((df$HEIGHT - mean(df$HEIGHT))^2)

  # ---- Mixed-effects -------------------------------------------------------
  if (!is.null(group_col)) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop(
        "Package 'nlme' is required for mixed-effects fitting. ",
        "Install it with: install.packages('nlme')"
      )
    }

    n_params   <- length(start)
    fixed_form <- if (n_params == 2) {
      stats::as.formula("a + b ~ 1")
    } else {
      stats::as.formula("a + b + c ~ 1")
    }

    random_form <- if (!is.null(nested_col)) {
      stats::as.formula("a ~ 1 | SITE_IDENTIFIER/unitreeid")
    } else {
      stats::as.formula("a ~ 1 | SITE_IDENTIFIER")
    }

    fit <- try(
      nlme::nlme(
        model   = formula_fit,
        data    = df,
        fixed   = fixed_form,
        start   = unlist(start),
        random  = random_form,
        control = nlme::nlmeControl(minScale = 1e-100, maxIter = 1000)
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      stop(
        "Mixed-effects model fitting failed for model '", model, "'. ",
        "Try a different model form or check your data."
      )
    }

    resid_df       <- data.frame(fit$residuals)
    r2_marg        <- 1 - sum(resid_df[["fixed"]]^2) / R2_bottom
    # Conditional R²: innermost residual column (matches FAIBCompiler r2_cond)
    inner_col      <- if (!is.null(nested_col)) "unitreeid" else "SITE_IDENTIFIER"
    r2_cond        <- 1 - sum(resid_df[[inner_col]]^2) / R2_bottom
    coefs          <- fit$coefficients$fixed

    # Extract site-level random effects (BLUPs for asymptote parameter a)
    re_site                  <- data.frame(fit$coefficients$random$SITE_IDENTIFIER)
    colnames(re_site)        <- sub("\\(Intercept\\)", "a", colnames(re_site))
    re_site$SITE_IDENTIFIER  <- rownames(re_site)
    rownames(re_site)        <- NULL

    out <- list(
      model                = model,
      fit                  = fit,
      coefficients         = coefs,
      r2_marg              = r2_marg,
      r2_cond              = r2_cond,
      random_effects_site  = re_site[, c("SITE_IDENTIFIER", "a"), drop = FALSE]
    )

    # Extract nested-group random effects when two-level nesting was used
    if (!is.null(nested_col)) {
      re_tree           <- data.frame(fit$coefficients$random$unitreeid)
      colnames(re_tree) <- sub("\\(Intercept\\)", "a", colnames(re_tree))
      re_tree$.key      <- rownames(re_tree)
      rownames(re_tree) <- NULL
      # nlme nested row names are formatted as "site_value/tree_value"
      re_tree$unitreeid <- sub(".*/", "", re_tree$.key)
      re_tree$.key      <- NULL
      out$random_effects_tree <- re_tree[, c("unitreeid", "a"), drop = FALSE]
    }

    return(out)

  # ---- Fixed effects -------------------------------------------------------
  } else {
    fit <- try(
      stats::nls(formula = formula_fit, data = df, start = start),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      stop(
        "Model fitting failed for model '", model, "'. ",
        "Try a different model form or check your data."
      )
    }

    r2_marg <- 1 - sum(stats::residuals(fit)^2) / R2_bottom
    coefs       <- stats::coef(fit)
  }

  list(
    model        = model,
    fit          = fit,
    coefficients = coefs,
    r2_marg      = r2_marg
  )
}


# Internal: apply one H-D equation given scalar (or per-row vectorised) coefs.
# Called row-by-row from ht_predict() via mapply so c may be NA for 2-param
# models and the switch will simply ignore it.
.hd_apply_coefs <- function(dbh, model, a, b, c = NA_real_) {
  switch(
    tolower(model),
    naslund  = 1.3 + dbh^2 / (a + b * dbh)^2,
    curtis   = 1.3 + a * (dbh / (1 + dbh))^b,
    logistic = 1.3 + a / (1 + b * exp(-c * dbh)),
    korf     = 1.3 + a * exp(-b * dbh^(-c)),
    weibull  = 1.3 + a * (1 - exp(-b * dbh^c)),
    richards = 1.3 + a * (1 - exp(-b * dbh))^c,
    {
      warning(".hd_apply_coefs: unknown model form '", model, "'.")
      NA_real_
    }
  )
}


#' Predict tree heights from a hierarchy of fitted H-D models
#'
#' @description
#' Applies a prioritised list of pre-fitted height-diameter coefficient tables
#' to a tree data frame, working from the most specific to the most general
#' model until every tree is matched.  Each tree records \emph{which level} of
#' the hierarchy supplied its model, making fallbacks explicit rather than
#' silent.
#'
#' @details
#' \strong{How the hierarchy works:}
#' \code{models} is a \emph{named} list of data frames, tried in order.  For
#' each level the function performs a left-join of the still-unmatched trees
#' against the coefficient table using whatever columns the coefficient table
#' contains (other than \code{model}, \code{a}, \code{b}, \code{c}).  Trees
#' that find a match at level \eqn{k} are predicted immediately and removed
#' from the unmatched pool; remaining trees cascade to level \eqn{k+1}.
#'
#' \strong{Coefficient table format:}
#' Each element of \code{models} must be a data frame with:
#' \itemize{
#'   \item One or more \emph{group columns} whose names match columns in
#'     \code{data} (e.g. \code{BEC_ZONE}, \code{SPECIES}, \code{SP0},
#'     \code{SP_TYPE}).  These are used for the join.  A level with
#'     \emph{no} group columns acts as a global catch-all and is applied to
#'     every remaining unmatched tree (useful as the last fallback).
#'   \item \code{model}: character -- equation form name (one of
#'     \code{"naslund"}, \code{"curtis"}, \code{"logistic"}, \code{"korf"},
#'     \code{"weibull"}, \code{"richards"}).
#'   \item \code{a}, \code{b}: numeric coefficient columns (required).
#'   \item \code{c}: numeric coefficient column (optional; omit for
#'     two-parameter forms such as \code{naslund} and \code{curtis}).
#' }
#' Each group-key combination must be unique within a table.
#'
#' \strong{Flag column:}
#' The \code{source_col} column in the returned data frame records the
#' \emph{name} of the list element that provided the model for each tree
#' (e.g. \code{"bec_species"}, \code{"sp0"}, \code{"sp_type"}).
#' Trees that could not be matched to any level receive \code{NA} in both the
#' height and source columns, and a warning is issued -- heights are never
#' silently imputed.
#'
#' @param data       data frame. Must contain the DBH column and all group
#'   columns referenced by any element of \code{models}.
#' @param models     named list of data frames, ordered from most specific to
#'   most general (see Details).
#' @param dbh_col    character. Name of the DBH column in \code{data}.
#'   Default \code{"DBH"}.
#' @param ht_col     character. Name of the predicted-height column added to
#'   the output. Default \code{"HT_PREDICTED"}.
#' @param source_col character. Name of the model-source flag column added to
#'   the output. Default \code{"HT_MODEL_SOURCE"}.
#'
#' @return \code{data} with two additional columns: \code{ht_col} (numeric,
#'   predicted heights in m) and \code{source_col} (character, name of the
#'   hierarchy level that supplied the model).  Unmatched trees have \code{NA}
#'   in both columns.
#'
#' @export
#' @seealso \code{\link{fit_hd_model}}, \code{\link{ht_from_dbh}}
#' @examples
#' set.seed(1)
#' trees <- data.frame(
#'   SPECIES = c("Pl", "Sx", "Hw", "Fd", "Bl"),
#'   SP0     = c("PL", "S",  "H",  "F",  "B"),
#'   DBH     = c(15,   22,   18,   30,   12)
#' )
#'
#' # Coefficient tables at two levels
#' coefs_sp <- data.frame(
#'   SPECIES = c("Pl", "Sx", "Hw"),
#'   model   = "naslund",
#'   a = c(1.2, 1.0, 1.5), b = c(0.04, 0.05, 0.03)
#' )
#' coefs_sp0 <- data.frame(
#'   SP0   = c("F", "B"),
#'   model = "naslund",
#'   a = c(1.1, 1.3), b = c(0.045, 0.035)
#' )
#'
#' result <- ht_predict(trees,
#'                      models = list(species = coefs_sp, sp0 = coefs_sp0))
#' result[, c("SPECIES", "DBH", "HT_PREDICTED", "HT_MODEL_SOURCE")]
ht_predict <- function(data,
                       models,
                       dbh_col    = "DBH",
                       ht_col     = "HT_PREDICTED",
                       source_col = "HT_MODEL_SOURCE") {

  if (!is.list(models) || is.null(names(models)) || any(names(models) == "")) {
    stop("'models' must be a named list (no unnamed elements).")
  }
  if (!dbh_col %in% names(data)) {
    stop("DBH column '", dbh_col, "' not found in data.")
  }

  out              <- as.data.frame(data, stringsAsFactors = FALSE)
  out[[ht_col]]    <- NA_real_
  out[[source_col]] <- NA_character_
  unmatched        <- rep(TRUE, nrow(out))   # TRUE = still needs a model

  reserved <- c("model", "a", "b", "c",
                 "n_obs", "r2_marg", "r2_cond", "fit_status")  # metadata from fit_hd_models_by_group()

  for (level_name in names(models)) {
    if (!any(unmatched)) break

    coef_df  <- as.data.frame(models[[level_name]], stringsAsFactors = FALSE)
    grp_cols <- setdiff(names(coef_df), reserved)

    # Validate: all group columns must exist in data
    missing_grp <- setdiff(grp_cols, names(out))
    if (length(missing_grp) > 0L) {
      warning("ht_predict: level '", level_name,
              "' references column(s) not found in data (",
              paste(missing_grp, collapse = ", "), ") -- skipping.")
      next
    }

    # Warn on duplicate group keys
    if (length(grp_cols) > 0L &&
        anyDuplicated(coef_df[, grp_cols, drop = FALSE]) > 0L) {
      warning("ht_predict: level '", level_name,
              "' has duplicate group-key rows -- only the first match per ",
              "tree will be used.")
    }

    # Subset unmatched trees; attach row index so it survives merge reordering
    idx          <- which(unmatched)
    sub          <- out[idx, , drop = FALSE]
    sub$.row_id  <- idx

    keep_coef <- intersect(c(grp_cols, reserved), names(coef_df))

    if (length(grp_cols) > 0L) {
      merged <- merge(sub, coef_df[, keep_coef, drop = FALSE],
                      by = grp_cols, all.x = TRUE, sort = FALSE)
    } else {
      # No group columns: apply the single catch-all row to every remaining tree
      if (nrow(coef_df) > 1L) {
        warning("ht_predict: level '", level_name,
                "' has no group columns but ", nrow(coef_df),
                " rows -- using only the first row.")
      }
      coef_row <- coef_df[1L, intersect(reserved, names(coef_df)), drop = FALSE]
      merged   <- cbind(sub,
                        coef_row[rep(1L, nrow(sub)), , drop = FALSE],
                        row.names = NULL)
    }

    has_match <- !is.na(merged$model) & !is.na(merged$a) & !is.na(merged$b)
    if (!any(has_match)) next

    m     <- merged[has_match, , drop = FALSE]
    c_vec <- if ("c" %in% names(m)) m$c else rep(NA_real_, nrow(m))

    ht_pred <- mapply(
      .hd_apply_coefs,
      dbh   = m[[dbh_col]],
      model = m$model,
      a     = m$a,
      b     = m$b,
      c     = c_vec,
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )

    row_ids                    <- m$.row_id
    out[[ht_col]][row_ids]     <- ht_pred
    out[[source_col]][row_ids] <- level_name
    unmatched[row_ids]         <- FALSE
  }

  n_unmatched <- sum(unmatched)
  if (n_unmatched > 0L) {
    warning("ht_predict: ", n_unmatched,
            " tree(s) could not be matched to any model level; ",
            "heights set to NA.")
  }

  out
}


#' Fit height-diameter models across SP0 groups with SP_TYPE fallback
#'
#' @description
#' Fits height-diameter allometric models to every SP0 species group found in
#' \code{data}, then automatically falls back to the broader SP_TYPE
#' (conifer/deciduous) level for any SP0 group that had too few trees or failed
#' to converge, following the two-level hierarchy used in the BC MoF PSP and
#' non-PSP compilation routines.
#'
#' For each group the function tries mixed-effects fitting first (when
#' \pkg{nlme} is available and enough distinct plots are present), then falls
#' back to fixed-effects (\code{nls}).  When \code{model} contains more than
#' one form, all are tried and the best by conditional \eqn{R^2} (mixed) or
#' marginal \eqn{R^2} (fixed) is retained.
#'
#' The returned \code{fits} list preserves the full model objects, including
#' \code{random_effects_site} (site-level BLUPs), so that predictions can
#' incorporate the plot-level random effect
#' (\eqn{a_{\text{total}} = a_{\text{fixed}} + u_{\text{site}}}).
#'
#' @param data        data frame of \strong{measured} trees (trees with both
#'   DBH and height recorded; broken-top trees must already be excluded).
#' @param sp0_col     character. Column name for the SP0 species group code.
#'   Default \code{"SPECIES_SP0"}.
#' @param sptype_col  character. Column name for the species type code
#'   (\code{"C"} = conifer, \code{"D"} = deciduous). Used only for the
#'   SP_TYPE fallback level.  Default \code{"SPECIES_SP_TYPE"}.
#' @param site_col    character. Column name for the site/plot identifier
#'   (\code{SITE_IDENTIFIER} in BC PSP data). Required for mixed-effects
#'   fitting.  Default \code{"SITE_IDENTIFIER"}.
#' @param tree_col    character. Column name for the individual tree identifier
#'   (\code{unitreeid} in BC PSP data). When present and trees appear more
#'   than once (multi-year data), a two-level nested random effect
#'   (\code{site / tree}) is used.  Default \code{"unitreeid"}.
#' @param height_col  character. Name of the height column. Default
#'   \code{"HEIGHT"}.
#' @param dbh_col     character. Name of the DBH column. Default \code{"DBH"}.
#' @param model       character vector. Model form(s) to attempt for each
#'   group.  When more than one is supplied, all are tried and the best is
#'   selected.  Default \code{c("naslund", "curtis")}.
#' @param min_n       integer. Minimum number of valid observations required
#'   to attempt fitting for an SP0 group.  Default \code{10L}, which suits
#'   small study-area datasets.  The BC MoF compilation routine
#'   (\code{DBH_Height_MEM.R}) uses \strong{1 000 observations} as its
#'   threshold before fitting a species- or SP0-level model from the
#'   provincial PSP + non-PSP pool; users with province-wide data may wish to
#'   set \code{min_n = 1000L} to match that standard.
#' @param min_plots   integer. Minimum number of distinct plots required to
#'   attempt mixed-effects fitting.  Default \code{5L}.  Increase this if
#'   you want tighter guarantees on the between-plot variance estimate (the
#'   BC MoF routine does not enforce an explicit plot minimum — it simply
#'   attempts the mixed-effects fit and falls back to fixed-effects on
#'   convergence failure).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{fits}}{Named list, one element per group (SP0 groups first,
#'       then SP_TYPE fallbacks). Each element is itself a list with:
#'       \itemize{
#'         \item \code{fit} — output of \code{\link{fit_hd_model}} (includes
#'           \code{coefficients}, \code{r2_marg}, \code{r2_cond},
#'           \code{random_effects_site}, etc.), or \code{NULL} if no model
#'           could be fitted.
#'         \item \code{method} — \code{"mixed-effects"},
#'           \code{"fixed-effects"}, or \code{"none"}.
#'         \item \code{level} — \code{"sp0"} or \code{"sp_type"}.
#'         \item \code{group} — the group identifier value (SP0 code or
#'           SP_TYPE code).
#'         \item \code{n} — number of training observations.
#'       }
#'     }
#'     \item{\code{summary}}{data frame with columns \code{Level},
#'       \code{Group}, \code{Method}, \code{n}, \code{Marg_R2}.}
#'   }
#'
#' @export
#' @seealso \code{\link{fit_hd_model}}, \code{\link{ht_from_dbh}}
#' @examples
#' \dontrun{
#' # Assumes psp_trees has been crosswalked to add SPECIES_SP0, SPECIES_SP_TYPE
#' library(CanWestAllometryCFS)
#' trees <- psp_trees
#' trees$SPECIES_CORR    <- species_correction(trees$SPECIES, trees$BEC_ZONE)
#' trees$SPECIES_SP0     <- bc_species_to_sp0(trees$SPECIES_CORR)
#' trees$SPECIES_SP_TYPE <- bc_species_to_sp_type(trees$SPECIES_CORR)
#'
#' measured  <- trees[!is.na(trees$HEIGHT) & !trees$BTOP, ]
#' hd_result <- fit_hd_models_by_group(measured)
#' print(hd_result$summary, row.names = FALSE)
#' }
fit_hd_models_by_group <- function(data,
                                    sp0_col    = "SPECIES_SP0",
                                    sptype_col = "SPECIES_SP_TYPE",
                                    site_col   = "SITE_IDENTIFIER",
                                    tree_col   = "unitreeid",
                                    height_col = "HEIGHT",
                                    dbh_col    = "DBH",
                                    model      = c("naslund", "curtis"),
                                    min_n      = 10L,
                                    min_plots  = 5L) {

  # --- Validate required columns -------------------------------------------
  required <- c(sp0_col, sptype_col, height_col, dbh_col)
  missing  <- setdiff(required, names(data))
  if (length(missing) > 0L)
    stop("Column(s) not found in data: ", paste(missing, collapse = ", "))

  has_site <- site_col %in% names(data)
  has_tree <- tree_col %in% names(data)

  # --- Filter valid rows -----------------------------------------------------
  keep <- !is.na(data[[dbh_col]])      &
          !is.na(data[[height_col]])   &
          data[[dbh_col]]    > 0       &
          data[[height_col]] > 1.3     &
          !is.na(data[[sp0_col]])      &
          !is.na(data[[sptype_col]])
  data_clean <- data[keep, , drop = FALSE]
  if (nrow(data_clean) == 0L)
    stop("No valid rows remain after filtering.")

  # --- Internal helper: fit one group (best model form, ME then FE) ----------
  .fit_one <- function(sub) {
    n_plots <- if (has_site) length(unique(sub[[site_col]])) else 0L
    best_fit <- NULL; best_method <- "none"; best_r2 <- -Inf

    # Try mixed-effects
    if (has_site && nrow(sub) >= min_n && n_plots >= min_plots &&
        requireNamespace("nlme", quietly = TRUE)) {
      has_rep <- has_tree && anyDuplicated(sub[[tree_col]]) > 0L
      for (mform in model) {
        fit <- tryCatch(
          fit_hd_model(sub, model = mform,
                       height_col = height_col, dbh_col = dbh_col,
                       group_col  = site_col,
                       nested_col = if (has_rep) tree_col else NULL),
          error = function(e) NULL)
        if (is.null(fit)) next
        r2 <- if (!is.null(fit$r2_cond)) fit$r2_cond else fit$r2_marg
        if (r2 > best_r2) { best_r2 <- r2; best_fit <- fit; best_method <- "mixed-effects" }
      }
    }
    # Fixed-effects fallback
    if (is.null(best_fit) && nrow(sub) >= min_n) {
      for (mform in model) {
        fit <- tryCatch(
          fit_hd_model(sub, model = mform,
                       height_col = height_col, dbh_col = dbh_col),
          error = function(e) NULL)
        if (is.null(fit)) next
        if (fit$r2_marg > best_r2) {
          best_r2 <- fit$r2_marg; best_fit <- fit; best_method <- "fixed-effects"
        }
      }
    }
    list(fit = best_fit, method = best_method, n = nrow(sub))
  }

  # --- Level 1: SP0 fits -----------------------------------------------------
  sp0_groups <- sort(unique(data_clean[[sp0_col]]))
  fits_sp0   <- lapply(sp0_groups, function(sp0) {
    sub <- data_clean[data_clean[[sp0_col]] == sp0, , drop = FALSE]
    c(.fit_one(sub), list(level = "sp0", group = sp0))
  })
  names(fits_sp0) <- sp0_groups

  # --- Level 2: SP_TYPE fallback for SP0 groups with no model ----------------
  sp0_no_model  <- vapply(fits_sp0, function(m) m$method == "none", logical(1))
  sptype_needed <- sort(unique(
    data_clean[[sptype_col]][
      data_clean[[sp0_col]] %in% names(sp0_no_model)[sp0_no_model]]
  ))
  fits_sptype <- lapply(sptype_needed, function(spt) {
    sub <- data_clean[data_clean[[sptype_col]] == spt, , drop = FALSE]
    c(.fit_one(sub), list(level = "sp_type", group = spt))
  })
  names(fits_sptype) <- sptype_needed

  # --- Combine and build summary ---------------------------------------------
  all_fits   <- c(fits_sp0, fits_sptype)
  summary_df <- do.call(rbind, lapply(all_fits, function(m) {
    data.frame(
      Level   = m$level,
      Group   = m$group,
      Method  = m$method,
      n       = m$n,
      Marg_R2 = if (!is.null(m$fit)) round(m$fit$r2_marg, 4L) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_df) <- NULL

  list(fits = all_fits, summary = summary_df)
}


#' Impute missing tree heights from fitted H-D models
#'
#' @description
#' Applies the hierarchical H-D model fits produced by
#' \code{\link{fit_hd_models_by_group}} to a tree data frame, filling missing
#' heights in two passes (SP0-level models first, then SP_TYPE fallbacks for
#' trees still unmatched).  Where a mixed-effects model was fitted, predictions
#' combine the fixed-effect coefficients with plot-level BLUPs
#' (\eqn{a_{\text{pred}} = a_{\text{fixed}} + u_{\text{site}}}).
#'
#' When \code{\link{fit_hd_models_by_group}} was called with a
#' \code{tree_col} (multi-year data), the fitted model contains both
#' plot-level and tree-level random effects.  Only the plot-level BLUP is
#' applied here — it localises the estimate to the plot's productivity level
#' and is available for every plot in the training data.  The tree-level BLUP
#' is not applied during imputation because imputation targets trees whose
#' height was never measured, so no tree-level BLUP exists for them.
#'
#' By default (\code{impute_btop = FALSE}), broken-top trees are excluded from
#' imputation: their height remains \code{NA} and they receive
#' \code{HT_FLAG = "btop"}.  Set \code{impute_btop = TRUE} to also estimate
#' total height for broken-top trees from DBH using the same H-D model.
#' This replicates BC MoF's compilation routine behaviour: all broken-top
#' trees receive a DBH-based height estimate, and every tree then receives an
#' \code{HT_PROJ} so biomass can be calculated in a single pass using
#' height-included equations rather than a DBH-only fallback.
#'
#' When \code{visit_col} is supplied (multi-year PSP data), a QC check
#' is applied: for each tree that has both measured and estimated heights
#' across visits, the function checks whether any estimated height is more
#' than 3 m higher than the measured height in the immediately following
#' visit.  Such trees are flagged \code{"ht_est_too_high"} in
#' \code{flag_col} because the discrepancy suggests an over-prediction.
#'
#' @param data       data frame. Must contain the columns named by
#'   \code{height_col}, \code{dbh_col}, \code{sp0_col}, \code{sptype_col},
#'   \code{site_col}, and \code{btop_col}.
#' @param hd_result  list returned by \code{\link{fit_hd_models_by_group}}.
#' @param height_col character. Name of the measured-height column.
#'   Default \code{"HEIGHT"}.
#' @param dbh_col    character. Name of the DBH column. Default \code{"DBH"}.
#' @param sp0_col    character. Name of the SP0 species group column.
#'   Default \code{"SPECIES_SP0"}.
#' @param sptype_col character. Name of the SP_TYPE column (\code{"C"}/\code{"D"}).
#'   Default \code{"SPECIES_SP_TYPE"}.
#' @param site_col   character. Name of the site/plot identifier column.
#'   Default \code{"SITE_IDENTIFIER"}.
#' @param btop_col   character. Name of the broken-top indicator column
#'   (\code{logical}). Default \code{"BTOP"}.
#' @param visit_col  character or \code{NULL}. Name of the visit/year column.
#'   When supplied, enables the multi-year QC check (see Details).
#'   Default \code{NULL}.
#' @param impute_btop logical. If \code{TRUE}, the H-D model is also used to
#'   estimate total height for broken-top trees (from DBH only, which is
#'   unaffected by the break). The resulting height represents what the tree
#'   would have reached had it remained intact, and is useful for biomass
#'   estimation. The \code{flag_col} column retains \code{"btop"} for these
#'   trees so they remain identifiable.  Default \code{FALSE}.
#' @param ht_col     character. Name of the output height column added to
#'   \code{data}. Default \code{"HT_PROJ"}.
#' @param source_col character. Name of the audit column recording the height
#'   source for each tree. Default \code{"HD_SOURCE"}.
#' @param flag_col   character. Name of the QC flag column added to
#'   \code{data}. Default \code{"HT_FLAG"}.
#'
#' @return \code{data} with three additional columns:
#'   \describe{
#'     \item{\code{ht_col}}{Numeric. Filled heights (m). Trees with a
#'       measured height retain that value; trees without a measured height
#'       receive a model-based prediction; broken-top trees remain \code{NA}
#'       unless \code{impute_btop = TRUE}.}
#'     \item{\code{source_col}}{Character. \code{"measured"} for trees with
#'       an observed height; \code{"naslund_mixed-effects_sp0"} etc. for
#'       model-imputed trees; \code{NA} for trees that could not be matched
#'       (no model for their SP0/SP_TYPE group, or broken top).}
#'     \item{\code{flag_col}}{Character. \code{"btop"} for broken-top trees;
#'       \code{"ht_est_too_high"} for trees where the estimated height exceeds
#'       the subsequent measured height by more than 3 m (multi-year data only);
#'       \code{NA} otherwise.}
#'   }
#'
#' @export
#' @seealso \code{\link{fit_hd_models_by_group}}, \code{\link{fit_hd_model}}
#' @examples
#' \dontrun{
#' library(CanWestAllometryCFS)
#' trees <- psp_trees
#' trees$SPECIES_CORR    <- species_correction(trees$SPECIES, trees$BEC_ZONE)
#' trees$SPECIES_SP0     <- bc_species_to_sp0(trees$SPECIES_CORR)
#' trees$SPECIES_SP_TYPE <- bc_species_to_sp_type(trees$SPECIES_CORR)
#'
#' measured  <- trees[!is.na(trees$HEIGHT) & !trees$BTOP, ]
#' hd_result <- fit_hd_models_by_group(measured)
#'
#' # Default: broken-top trees left with HT_PROJ = NA
#' trees <- ht_impute(trees, hd_result)
#' table(trees$HT_FLAG, useNA = "ifany")
#'
#' # BC MoF compilation routine: estimate height for broken-top trees too.
#' # Setting impute_btop = TRUE means every tree receives an HT_PROJ,
#' # allowing single-pass biomass estimation with height-included equations.
#' trees <- ht_impute(trees, hd_result, impute_btop = TRUE)
#' trees$BIOMASS_KG <- biomass_tree(trees$SPECIES_NAME, trees$DBH,
#'                                   height = trees$HT_PROJ)
#'
#' # Multi-year data: also pass visit_col to enable the QC check
#' # trees <- ht_impute(trees, hd_result, impute_btop = TRUE,
#' #                    visit_col = "VISIT_NUMBER")
#' }
ht_impute <- function(data,
                      hd_result,
                      height_col  = "HEIGHT",
                      dbh_col     = "DBH",
                      sp0_col     = "SPECIES_SP0",
                      sptype_col  = "SPECIES_SP_TYPE",
                      site_col    = "SITE_IDENTIFIER",
                      btop_col    = "BTOP",
                      visit_col   = NULL,
                      impute_btop = FALSE,
                      ht_col      = "HT_PROJ",
                      source_col  = "HD_SOURCE",
                      flag_col    = "HT_FLAG") {

  # --- Validate inputs -------------------------------------------------------
  required <- c(height_col, dbh_col, sp0_col, sptype_col, btop_col)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L)
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  if (!is.null(visit_col) && !visit_col %in% names(data))
    stop("visit_col '", visit_col, "' not found in data.")
  if (!is.list(hd_result) || !all(c("fits", "summary") %in% names(hd_result)))
    stop("'hd_result' must be the list returned by fit_hd_models_by_group().")

  has_site <- site_col %in% names(data)

  # --- Initialise output columns ---------------------------------------------
  data[[ht_col]]     <- data[[height_col]]
  data[[source_col]] <- ifelse(!is.na(data[[height_col]]), "measured",
                               NA_character_)
  data[[flag_col]]   <- NA_character_

  # Flag broken-top trees immediately
  data[[flag_col]][data[[btop_col]]] <- "btop"

  # --- Internal BLUP lookup --------------------------------------------------
  .blup_a <- function(site_ids, re_site) {
    if (is.null(re_site) || !has_site)
      return(rep(0, length(site_ids)))
    b <- re_site$a[match(site_ids, re_site$SITE_IDENTIFIER)]
    ifelse(is.na(b), 0, b)
  }

  # --- Imputation loop (SP0 models run before SP_TYPE fallbacks) -------------
  for (m in hd_result$fits) {
    if (is.null(m$fit)) next
    cf      <- m$fit$coefficients
    re_site <- m$fit$random_effects_site
    model   <- m$fit$model

    if (m$level == "sp0") {
      # Normal pass: missing height AND not broken top
      missing <- is.na(data[[height_col]]) & !data[[btop_col]] &
                 !is.na(data[[sp0_col]]) & data[[sp0_col]] == m$group
      # Optional broken-top pass: broken top, not yet imputed
      if (impute_btop) {
        missing_btop <- data[[btop_col]] & is.na(data[[ht_col]]) &
                        !is.na(data[[sp0_col]]) & data[[sp0_col]] == m$group
        missing <- missing | missing_btop
      }
    } else {
      missing <- is.na(data[[ht_col]]) & !data[[btop_col]] &
                 !is.na(data[[sptype_col]]) & data[[sptype_col]] == m$group
      if (impute_btop) {
        missing_btop <- data[[btop_col]] & is.na(data[[ht_col]]) &
                        !is.na(data[[sptype_col]]) & data[[sptype_col]] == m$group
        missing <- missing | missing_btop
      }
    }
    if (!any(missing)) next

    site_ids <- if (has_site) data[[site_col]][missing] else
                  rep(NA_character_, sum(missing))
    a_pred   <- cf[["a"]] + .blup_a(site_ids, re_site)
    b_val    <- cf[["b"]]
    c_val    <- if ("c" %in% names(cf)) cf[["c"]] else NULL

    data[[ht_col]][missing] <- round(
      ht_from_dbh(data[[dbh_col]][missing], model = model,
                  a = a_pred, b = b_val, c = c_val),
      1L)
    data[[source_col]][missing] <- paste0(model, "_", m$method, "_", m$level)
  }

  # --- Multi-year QC: flag estimated heights that exceed next measured -------
  # Mirrors FAIBCompiler pspHT(): for trees with both measured and estimated
  # heights across visits, flag rows where estimated > next measured by > 3 m.
  if (!is.null(visit_col) && has_site && "unitreeid" %in% names(data)) {
    tree_id_col <- "unitreeid"
    # Sort by tree and visit
    ord <- order(data[[tree_id_col]], data[[visit_col]])
    d   <- data[ord, ]

    # For each tree, find rows that follow a source change (est → measured)
    src   <- d[[source_col]]
    ht    <- d[[ht_col]]
    tid   <- d[[tree_id_col]]

    n <- nrow(d)
    flag_idx <- integer(0)
    for (i in seq_len(n - 1L)) {
      if (tid[i] == tid[i + 1L] &&
          src[i] != "measured" && !is.na(src[i]) &&
          src[i + 1L] == "measured" && !is.na(ht[i]) && !is.na(ht[i + 1L]) &&
          (ht[i] - ht[i + 1L]) > 3) {
        flag_idx <- c(flag_idx, i)
      }
    }
    if (length(flag_idx) > 0L) {
      orig_rows <- ord[flag_idx]
      data[[flag_col]][orig_rows] <- "ht_est_too_high"
    }
  }

  data
}


# Internal: partial derivative dH/da for each model form --------------------
# Used by calibrate_hd_blup() to build the Z vector for the BLUP calculation.
# The random effect in fit_hd_model() is always placed on parameter 'a', so
# the gradient w.r.t. the random effect equals dH/da evaluated at u = 0.

.dH_da <- function(model, dbh, a, b, c = NULL) {
  switch(model,
    naslund  = -2 * dbh^2 / (a + b * dbh)^3,
    curtis   = (dbh / (1 + dbh))^b,
    logistic = 1 / (1 + b * exp(-c * dbh)),
    korf     = exp(-b * dbh^(-c)),
    weibull  = 1 - exp(-b * dbh^c),
    richards = (1 - exp(-b * dbh))^c,
    stop("Unknown model form: ", model)
  )
}


#' BLUP calibration of a mixed-effects H-D model to a new plot
#'
#' @description
#' Given a mixed-effects height-diameter model fitted with
#' \code{\link{fit_hd_model}} and one or more trees with \emph{measured}
#' height and DBH from a plot that was \strong{not} in the training data,
#' estimates the plot-level random effect using empirical Best Linear Unbiased
#' Prediction (BLUP) and returns a calibrated coefficient set for that plot.
#'
#' The calibrated \code{a} can then be passed to \code{\link{ht_from_dbh}}
#' to predict heights for the remaining unmeasured trees in the same plot with
#' a plot-specific curve rather than the population average.
#'
#' \strong{When to use this function:}
#' \itemize{
#'   \item You have a \emph{fixed, frozen} reference model (e.g. a published
#'     regional coefficient table) and new plots arrive that were not in the
#'     original training data. You do not have access to the training data and
#'     therefore cannot refit.
#'   \item You are working with a single-stand cruise (not a PSP remeasurement
#'     database) and want to apply a regional model with local calibration.
#'   \item The model is updated on a long cycle and you need to handle new
#'     plots between updates without triggering a full refit.
#' }
#'
#' \strong{When NOT to use this function:}
#' If you have access to the training data and new measurements have arrived,
#' the correct approach is to pool all data and refit the model with
#' \code{\link{fit_hd_model}}. The resulting mixed-effects fit will
#' automatically embed BLUPs for every plot that contributes measured trees,
#' and \code{predict(fit$fit, level = 1)} will give calibrated predictions
#' without any additional step. \code{calibrate_hd_blup} is redundant in
#' that workflow.
#'
#' \strong{Why this matters:} The mixed-effects model captures
#' plot-to-plot variation in site productivity via a random effect on the
#' asymptote parameter \eqn{a}. For plots not in the training data only the
#' fixed-effects (population-average) curve is available. Even one
#' measured height-DBH pair from the new plot allows back-calculating an
#' estimate of that plot's random effect, substantially reducing prediction
#' error on productive or unproductive sites.
#'
#' \strong{Statistical method:} A first-order linearisation of the nonlinear
#' model around the population mean gives the Henderson BLUP equation
#' (scalar random effect form):
#'
#' \deqn{\hat{u} = \frac{\sigma^2_u \sum_i z_i r_i}
#'                       {\sigma^2_\varepsilon + \sigma^2_u \sum_i z_i^2}}
#'
#' where \eqn{z_i = \partial H / \partial a} evaluated at the fixed-effect
#' coefficients, \eqn{r_i} is the residual from the population-average
#' prediction, \eqn{\sigma^2_u} is the random-effect variance, and
#' \eqn{\sigma^2_\varepsilon} is the residual variance -- both taken from the
#' fitted \code{nlme} object.
#'
#' The calibrated asymptote is \eqn{\hat{a}_{\text{plot}} = a + \hat{u}}.
#'
#' @param fit        A list returned by \code{\link{fit_hd_model}} with
#'   \code{group_col} supplied (i.e. a mixed-effects fit). The list must
#'   contain elements \code{fit} (an \code{nlme} object), \code{coefficients},
#'   and \code{model}.
#' @param calib_data A data frame of trees from the \strong{new plot} that
#'   have both DBH and height measured. Must contain the columns named by
#'   \code{dbh_col} and \code{height_col}.
#' @param dbh_col    character. Name of the DBH column in \code{calib_data}.
#'   Default \code{"DBH"}.
#' @param height_col character. Name of the height column in \code{calib_data}.
#'   Default \code{"HEIGHT"}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{coefficients}}{Named numeric vector with calibrated
#'       coefficients for the target plot. The \code{a} element equals
#'       \eqn{a_{\text{fixed}} + \hat{u}}; all other coefficients are
#'       unchanged. Pass directly to \code{\link{ht_from_dbh}}.}
#'     \item{\code{u_hat}}{Numeric. Estimated random effect \eqn{\hat{u}}.
#'       Positive values indicate a more productive site than average;
#'       negative values indicate a less productive site.}
#'     \item{\code{model}}{character. Equation form (inherited from
#'       \code{fit}).}
#'     \item{\code{n_calib}}{Integer. Number of calibration trees used.}
#'   }
#'
#' @export
#' @seealso \code{\link{fit_hd_model}}, \code{\link{ht_from_dbh}}
#' @references
#'   Calama, R. and Montero, G. (2004). Interregional nonlinear height-diameter
#'   model with random coefficients for stone pine in Spain.
#'   \emph{Canadian Journal of Forest Research}, 34, 150--163.
#'
#'   Vonesh, E.F. and Chinchilli, V.M. (1997).
#'   \emph{Linear and Nonlinear Models for the Analysis of Repeated Measurements}.
#'   Marcel Dekker, New York.
#' @examples
#' \dontrun{
#' # 1. Fit a mixed-effects model on training data
#' set.seed(42)
#' n_plots  <- 12
#' plot_ids <- rep(paste0("P", 1:n_plots), each = 25)
#' dbh      <- runif(n_plots * 25, 5, 60)
#' a_plot   <- rep(rnorm(n_plots, mean = 1.5, sd = 0.3), each = 25)
#' height   <- ht_naslund(dbh, a = a_plot, b = 0.05) + rnorm(n_plots * 25, 0, 1)
#' train    <- data.frame(PLOT_ID = plot_ids, DBH = dbh, HEIGHT = height)
#'
#' fit <- fit_hd_model(train, model = "naslund", group_col = "PLOT_ID")
#'
#' # 2. New plot: a few measured trees (the calibration set)
#' calib <- data.frame(
#'   DBH    = c(15, 28, 42),
#'   HEIGHT = ht_naslund(c(15, 28, 42), a = 1.9, b = 0.05) + rnorm(3, 0, 0.5)
#' )
#'
#' # 3. Calibrate to the new plot
#' calib_result <- calibrate_hd_blup(fit, calib)
#' cat("Population-average a:", fit$coefficients["a"], "\n")
#' cat("Calibrated a        :", calib_result$coefficients["a"], "\n")
#' cat("Estimated u_hat     :", calib_result$u_hat, "\n")
#'
#' # 4. Predict unmeasured trees using the calibrated curve
#' new_dbh <- c(10, 20, 35, 50)
#' cf      <- calib_result$coefficients
#' ht_from_dbh(new_dbh, model = "naslund", a = cf["a"], b = cf["b"])
#' }
calibrate_hd_blup <- function(fit,
                               calib_data,
                               dbh_col    = "DBH",
                               height_col = "HEIGHT") {

  # --- Validate inputs -------------------------------------------------------
  if (!is.list(fit) || !all(c("fit", "coefficients", "model") %in% names(fit))) {
    stop("'fit' must be a list returned by fit_hd_model().")
  }
  if (!inherits(fit$fit, "nlme")) {
    stop(
      "BLUP calibration requires a mixed-effects model. ",
      "Refit with fit_hd_model(..., group_col = <plot column>)."
    )
  }

  missing_cols <- setdiff(c(dbh_col, height_col), names(calib_data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in calib_data: ",
         paste(missing_cols, collapse = ", "))
  }

  dbh_c <- calib_data[[dbh_col]]
  ht_c  <- calib_data[[height_col]]

  keep <- !is.na(dbh_c) & !is.na(ht_c) & dbh_c > 0 & ht_c > 1.3
  if (sum(keep) < 1L) {
    stop("No valid calibration trees (need DBH > 0 and HEIGHT > 1.3 m).")
  }
  dbh_c <- dbh_c[keep]
  ht_c  <- ht_c[keep]

  # --- Extract fixed-effect coefficients and variance components -------------
  cf      <- fit$coefficients
  a_fixed <- cf[["a"]]
  b_fixed <- cf[["b"]]
  c_fixed <- if ("c" %in% names(cf)) cf[["c"]] else NULL
  model   <- fit$model

  # Variance components from the nlme object
  vc         <- nlme::VarCorr(fit$fit)
  sigma_u_sq <- as.numeric(vc["a", "Variance"])
  sigma_e_sq <- as.numeric(vc["Residual", "Variance"])

  # --- BLUP formula ----------------------------------------------------------
  # Linearise around the population mean (u = 0):
  #   z_i = dH/da evaluated at fixed-effect coefficients
  z <- .dH_da(model, dbh_c, a_fixed, b_fixed, c_fixed)

  # Residuals from population-average prediction
  r <- ht_c - ht_from_dbh(dbh_c, model, a = a_fixed, b = b_fixed, c = c_fixed)

  # Scalar BLUP (Henderson equation for a single random effect):
  #   u_hat = sigma_u^2 * sum(z * r) / (sigma_e^2 + sigma_u^2 * sum(z^2))
  u_hat <- (sigma_u_sq * sum(z * r)) / (sigma_e_sq + sigma_u_sq * sum(z^2))

  # --- Return calibrated coefficients ----------------------------------------
  cf_calib    <- cf
  cf_calib[["a"]] <- a_fixed + u_hat

  list(
    coefficients = cf_calib,
    u_hat        = u_hat,
    model        = model,
    n_calib      = sum(keep)
  )
}
