# Portions of this file are adapted from the FAIBBase R package
# (https://github.com/bcgov/FAIBBase), file R/biomassCalculator.R.
# Copyright 2019 Province of British Columbia. Original author: Yong Luo.
# Licensed under the Apache License, Version 2.0; see the LICENSE file.
# Modified 2026 by Natural Resources Canada: restructured into vectorised
# biomass_tree() / biomass_components(), added component-wise output,
# height-optional equation selection, and roxygen documentation.

# Aboveground tree biomass using Lambert et al. (2005) and Ung et al. (2008)
# power allometric equations.  Exported functions: biomass_tree(),
# biomass_components().

# Internal lookup helper ------------------------------------------------------

.biomass_lookup <- function(species, paper_source, height_included) {
  sp_lower <- tolower(trimws(species))
  coefs    <- biomass_coefs[
    biomass_coefs$species         == sp_lower      &
    biomass_coefs$paper_source    == paper_source  &
    biomass_coefs$height_included == height_included,
  ]

  # If Ung2008 not available for this species, fall back to Lambert2005
  if (nrow(coefs) == 0 && paper_source == "Ung2008") {
    coefs <- biomass_coefs[
      biomass_coefs$species         == sp_lower       &
      biomass_coefs$paper_source    == "Lambert2005"  &
      biomass_coefs$height_included == height_included,
    ]
  }

  # For height-included path only: fall back to "generic" coefficients
  if (nrow(coefs) == 0 && height_included) {
    ps_try <- paper_source
    coefs  <- biomass_coefs[
      biomass_coefs$species         == "generic"  &
      biomass_coefs$paper_source    == ps_try      &
      biomass_coefs$height_included == TRUE,
    ]
    if (nrow(coefs) == 0) {
      coefs <- biomass_coefs[
        biomass_coefs$species         == "generic"      &
        biomass_coefs$paper_source    == "Lambert2005"  &
        biomass_coefs$height_included == TRUE,
      ]
    }
    attr(coefs, "used_generic") <- TRUE
  }

  coefs
}


# Core calculation (single tree, returns named numeric vector) ----------------

.biomass_calc <- function(dbh, height, coefs) {
  if (nrow(coefs) == 0 || is.na(dbh) || dbh <= 0) {
    return(c(wood = NA_real_, bark = NA_real_,
             branches = NA_real_, foliage = NA_real_))
  }
  height_included <- coefs$height_included[1]
  out <- numeric(4)
  names(out) <- c("wood", "bark", "branches", "foliage")
  for (comp in names(out)) {
    row <- coefs[coefs$component == comp, ]
    if (nrow(row) == 0) { out[comp] <- NA_real_; next }
    if (height_included) {
      a3 <- row$a3
      ht <- if (is.na(a3) || a3 == 0) 1 else height^a3
      out[comp] <- row$a1 * dbh^row$a2 * ht
    } else {
      out[comp] <- row$a1 * dbh^row$a2
    }
  }
  out
}


# Exported functions ----------------------------------------------------------

#' Compute total aboveground tree biomass
#'
#' @description
#' Computes total aboveground dry biomass (kg) for one or more trees using
#' the power allometric equations of Lambert et al. (2005) or Ung et al.
#' (2008). Biomass is the sum of wood, bark, branches, and foliage components.
#'
#' @details
#' \strong{Equations:}
#' \itemize{
#'   \item DBH-only (\code{height = NULL} or \code{height[i]} is \code{NA}):
#'     \eqn{B_{component} = a_1 \cdot DBH^{a_2}}
#'   \item DBH + height (\code{height[i]} is not \code{NA}):
#'     \eqn{B_{component} = a_1 \cdot DBH^{a_2} \cdot H^{a_3}}
#' }
#' Total biomass is summed across all four components.
#' When \code{height} is supplied as a vector, the equation is selected
#' per tree: trees with a non-\code{NA} height use the DBH + height
#' equations; trees with \code{NA} height automatically fall back to the
#' DBH-only equations.  This means a mixed dataset (some trees with heights,
#' some without) can be passed in a single call.
#'
#' \strong{Species names} must be lowercase common names matching those in
#' \code{\link{biomass_coefs}} (e.g. \code{"lodgepole pine"},
#' \code{"western redcedar"}). Use \code{sort(unique(biomass_coefs$species))}
#' to list all recognised names.
#'
#' \strong{Ung 2008 availability:} Only a subset of species have distinct
#' Ung 2008 coefficients. When \code{paper_source = "Ung2008"} is requested
#' for a species with Lambert 2005 coefficients only, Lambert 2005 values are
#' used silently.
#'
#' \strong{Unknown species:} For unrecognised species with
#' \code{height = NULL}, \code{NA} is returned with a warning. For the
#' height-included path, the generic softwood/hardwood average coefficients
#' from the source paper are used with a warning.
#'
#' @param species    character vector. Common species names (case-insensitive).
#' @param dbh        numeric vector. Diameter at breast height (cm). Must be
#'   the same length as \code{species}.
#' @param height     numeric vector or \code{NULL}. Total tree height (m). If
#'   supplied as a vector, trees where \code{height[i]} is not \code{NA}
#'   use the DBH + height equations; trees where \code{height[i]} is
#'   \code{NA} automatically fall back to the DBH-only equations. Pass
#'   \code{NULL} (default) to use DBH-only equations for all trees.
#'   Must be the same length as \code{species} if not \code{NULL}.
#' @param paper_source character. One of \code{"Lambert2005"} (default) or
#'   \code{"Ung2008"}.
#'
#' @return numeric vector of total aboveground dry biomass (kg), same length
#'   as \code{dbh}. \code{NA} is returned for trees where no coefficients
#'   can be found (DBH-only path only; see Details).
#'
#' @export
#' @seealso \code{\link{biomass_components}}, \code{\link{biomass_coefs}}
#' @references
#'   Lambert M-C, Ung C-H, Raulier F (2005). Canadian national tree aboveground
#'   biomass equations. \emph{Canadian Journal of Forest Research}, 35(8),
#'   1996–2018.
#'
#'   Ung C-H, Bernier P, Guo X-J (2008). Canadian national biomass equations:
#'   new parameter estimates that include British Columbia data.
#'   \emph{Canadian Journal of Forest Research}, 38(5), 1123–1132.
#'
#' @examples
#' # DBH-only (common PSP scenario without measured heights)
#' biomass_tree(species = c("lodgepole pine", "white spruce", "trembling aspen"),
#'              dbh     = c(15, 22, 18))
#'
#' # DBH + height (more accurate)
#' biomass_tree(species = c("lodgepole pine", "white spruce"),
#'              dbh     = c(15, 22),
#'              height  = c(12, 18))
#'
#' # Compare Lambert2005 vs Ung2008 for lodgepole pine
#' biomass_tree("lodgepole pine", 20, paper_source = "Lambert2005")
#' biomass_tree("lodgepole pine", 20, paper_source = "Ung2008")
biomass_tree <- function(species, dbh, height = NULL,
                         paper_source = "Lambert2005") {
  paper_source  <- match.arg(paper_source, c("Lambert2005", "Ung2008"))
  height_supplied <- !is.null(height)
  n             <- length(dbh)

  if (length(species) == 1L) species <- rep(species, n)
  if (length(species) != n) {
    stop("'species' and 'dbh' must be the same length (or 'species' scalar).")
  }
  if (height_supplied && length(height) != n) {
    stop("'height' must be the same length as 'dbh' when supplied.")
  }

  result <- numeric(n)

  for (i in seq_len(n)) {
    sp         <- tolower(trimws(species[i]))
    # Use height-included equations only when height[i] is non-NA
    height_incl <- height_supplied && !is.na(height[i])
    coefs <- .biomass_lookup(sp, paper_source, height_incl)

    if (nrow(coefs) == 0) {
      warning("No coefficients found for species '", species[i],
              "'. Returning NA.")
      result[i] <- NA_real_
      next
    }
    if (isTRUE(attr(coefs, "used_generic"))) {
      warning("Species '", species[i],
              "' not recognised; using generic coefficients.")
    }

    h        <- if (height_incl) height[i] else NA_real_
    comp_b   <- .biomass_calc(dbh[i], h, coefs)
    result[i] <- sum(comp_b, na.rm = FALSE)
  }

  result
}


#' Compute aboveground tree biomass by component
#'
#' @description
#' Returns a data frame of dry biomass (kg) broken down into wood, bark,
#' branches, and foliage for each tree. Uses the same equations and argument
#' conventions as \code{\link{biomass_tree}}.
#'
#' @inheritParams biomass_tree
#'
#' @return A data frame with \code{length(dbh)} rows and five columns:
#'   \code{wood}, \code{bark}, \code{branches}, \code{foliage}, and
#'   \code{total} (sum of all four components). All values are dry biomass in
#'   kilograms.
#'
#' @export
#' @seealso \code{\link{biomass_tree}}, \code{\link{biomass_coefs}}
#' @examples
#' biomass_components(
#'   species = c("douglas-fir", "western hemlock", "western redcedar"),
#'   dbh     = c(30, 25, 40),
#'   height  = c(28, 22, 35)
#' )
biomass_components <- function(species, dbh, height = NULL,
                               paper_source = "Lambert2005") {
  paper_source    <- match.arg(paper_source, c("Lambert2005", "Ung2008"))
  height_supplied <- !is.null(height)
  n               <- length(dbh)

  if (length(species) == 1L) species <- rep(species, n)
  if (length(species) != n) {
    stop("'species' and 'dbh' must be the same length (or 'species' scalar).")
  }
  if (height_supplied && length(height) != n) {
    stop("'height' must be the same length as 'dbh' when supplied.")
  }

  out <- data.frame(wood     = numeric(n),
                    bark     = numeric(n),
                    branches = numeric(n),
                    foliage  = numeric(n),
                    total    = numeric(n))

  for (i in seq_len(n)) {
    sp          <- tolower(trimws(species[i]))
    height_incl <- height_supplied && !is.na(height[i])
    coefs <- .biomass_lookup(sp, paper_source, height_incl)

    if (nrow(coefs) == 0) {
      warning("No coefficients found for species '", species[i],
              "'. Returning NA.")
      out[i, ] <- NA_real_
      next
    }
    if (isTRUE(attr(coefs, "used_generic"))) {
      warning("Species '", species[i],
              "' not recognised; using generic coefficients.")
    }

    h      <- if (height_incl) height[i] else NA_real_
    comp_b <- .biomass_calc(dbh[i], h, coefs)

    out$wood[i]     <- comp_b["wood"]
    out$bark[i]     <- comp_b["bark"]
    out$branches[i] <- comp_b["branches"]
    out$foliage[i]  <- comp_b["foliage"]
    out$total[i]    <- sum(comp_b)
  }

  out
}
