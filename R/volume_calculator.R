# Portions of this file are adapted from the FAIBBase R package
# (https://github.com/bcgov/FAIBBase), files R/treeVolCalculator.R,
# R/DIB_ICalculator.R and R/treeProfile.R.
# Copyright 2019 Province of British Columbia. Original author: Yong Luo.
# Licensed under the Apache License, Version 2.0; see the LICENSE file.
# Modified 2026 by Natural Resources Canada: taper coefficients moved out of
# the code into package datasets, KBEC/KFIZ3 forms exposed through a single
# tree_volume() interface, log-section and stem-profile helpers added, and
# roxygen documentation written.

# Kozak variable-exponent taper equations for BC tree species.
# Two stratification systems are supported:
#   KBEC  -- Kozak (2002 BC MoF internal / 2004 published) with BEC-zone
#             coefficients for 16 species x 13 BEC zones.
#   KFIZ3 -- earlier Kozak sqrt-transform form with Forest Inventory Zone
#             (FIZ) coefficients for 16 species x 12 FIZ zones.
# Exported functions: tree_volume(), tree_volume_section(), tree_profile().

# Internal Kozak (2002/2004) KBEC taper equation ----------------------------

#' @keywords internal
.kozak_kbec_dib <- function(height_i, dbh, height, coefs) {
  # Returns inside-bark diameter (cm) at height_i (m) for a single tree.
  # coefs: named numeric vector with B1:B9 and ERR.
  # Returns 0 at the tree top; never negative.
  zi  <- height_i / height
  p   <- 1.3 / height
  qi  <- 1 - zi^(1/3)
  xi  <- qi / (1 - p^(1/3))
  dh  <- dbh / height

  ex <- coefs[["B4"]] * zi^4     +
        coefs[["B5"]] / exp(dh)  +
        coefs[["B6"]] * xi^0.1   +
        coefs[["B7"]] / dbh      +
        coefs[["B8"]] * height^qi +
        coefs[["B9"]] * xi

  dib <- coefs[["ERR"]] *
         coefs[["B1"]] * (dbh^coefs[["B2"]]) *
         (height^coefs[["B3"]]) *
         (xi^ex)

  # At the tree top xi = 0; guard against NaN from 0^0 or 0^negative
  dib[!is.finite(dib) | xi <= 0] <- 0
  dib[dib < 0] <- 0
  dib
}


# Internal Kozak KFIZ3 (sqrt-transform) taper equation ----------------------

#' @keywords internal
.kozak_kfiz3_dib <- function(height_i, dbh, height, coefs) {
  # Returns inside-bark diameter (cm) at height_i (m) for a single tree.
  # coefs: named numeric vector with P, A0:A7.
  # KFIZ3 equation (Kozak sqrt-transform form):
  #   X = (1 - sqrt(Z)) / (1 - sqrt(P))  where P is a fitted inflection
  #   EX = A3*Z^2 + A4*log(Z+0.001) + A5*sqrt(X) + A6*(DBH/H) + A7*exp(Z)
  #   DIB = A0 * DBH^A1 * A2^DBH * X^EX
  z  <- height_i / height
  p  <- coefs[["P"]]
  x  <- (1 - sqrt(z)) / (1 - sqrt(p))

  ex <- coefs[["A3"]] * z^2 +
        coefs[["A4"]] * log(z + 0.001) +
        coefs[["A5"]] * sqrt(pmax(x, 0)) +
        coefs[["A6"]] * (dbh / height) +
        coefs[["A7"]] * exp(z)

  dib <- coefs[["A0"]] * (dbh^coefs[["A1"]]) * (coefs[["A2"]]^dbh) * (x^ex)

  # At tree top Z = 1, X = 0 -> DIB = 0
  dib[!is.finite(dib) | x <= 0] <- 0
  dib[dib < 0] <- 0
  dib
}


# Internal per-tree volume calculator ----------------------------------------
# dib_fn: a closure function(height_i) -> dib, already bound to the tree's
# dbh, height, and coefficient table row by the calling public function.

#' @keywords internal
.tree_vol_single <- function(dbh, height, volume_type, stump_height,
                             utop_dib, btop_height, dib_fn) {
  vcons <- pi / 40000  # converts DIB (cm) and height (m) -> m3

  # DIB at stump height (forced >= DIB at breast height for short stumps)
  dib_stump <- dib_fn(stump_height)
  dib_bh    <- dib_fn(1.3)
  if (dib_bh > dib_stump) dib_stump <- dib_bh

  vol_stump <- vcons * stump_height * dib_stump^2  # cylinder

  if (volume_type == "STUMP") return(vol_stump)

  # Respect broken top for all stem volume types: limit integration to upper
  upper <- if (!is.na(btop_height) && btop_height < height) btop_height else height

  # Integrate 10 cm slices from stump_height to upper using Smalian's formula
  heights <- seq(stump_height, upper, by = 0.1)
  if (tail(heights, 1) < upper) heights <- c(heights, upper)

  dibs <- dib_fn(heights)
  # Apply breast-height floor to all heights below 1.3 m, matching FAIBBase:
  # treeprofiledata[HT_I < breastHeight & DIB_I < DIB_BH, DIB_I := DIB_BH]
  below_bh <- heights < 1.3
  if (any(below_bh)) dibs[below_bh] <- pmax(dibs[below_bh], dib_bh)
  dibs[!is.finite(dibs) | dibs < 0] <- 0

  # Smalian: average of cross-section areas at slice ends x slice length
  areas      <- vcons * dibs^2
  lengths    <- diff(heights)
  vol_slices <- (head(areas, -1) + tail(areas, -1)) / 2 * lengths
  vol_stem   <- sum(vol_slices)  # above-stump volume up to upper

  if (volume_type == "WSV") return(vol_stump + vol_stem)

  # MER / NMR: check each slice by its top-of-slice DIB (matches FAIBBase
  # DIB_I_next >= UTOPDIB criterion; handles non-monotone taper correctly)
  slice_dib_tops <- dibs[-1]
  merch          <- slice_dib_tops >= utop_dib

  if (volume_type == "MER") return(max(0, sum(vol_slices[merch])))
  max(0, sum(vol_slices[!merch]))  # NMR
}


# Public function -------------------------------------------------------------

#' Calculate individual tree volume using Kozak taper equations
#'
#' @description
#' Computes individual tree volume (whole-stem, merchantable, non-merchantable,
#' or stump) using a Kozak variable-exponent taper equation.  Two taper
#' equation systems are supported via the \code{taper_eq} argument:
#'
#' \describe{
#'   \item{\code{"KBEC"} (default)}{Kozak's (2002 BC MoF internal report;
#'     published as Kozak 2004) variable-exponent taper equation with
#'     coefficients stratified by BC Biogeoclimatic Ecosystem Classification
#'     (BEC) zone and species.  Covers 16 species × 13 BEC zones.
#'     \strong{Use for any modern BC inventory data (VRI, current PSP data)
#'     where BEC zone is known.}  This is the preferred equation for
#'     post-2000 BC inventory work and matches the default in
#'     FAIBBase/FAIBCompiler.}
#'   \item{\code{"KFIZ3"}}{Kozak variable-exponent taper equation using a
#'     sqrt-transform shape variable (\eqn{X = (1-\sqrt{Z})/(1-\sqrt{P})})
#'     with coefficients stratified by the older BC Forest Inventory Zone
#'     (FIZ) system and species.  Covers 16 species × 12 FIZ zones (A–L),
#'     grouped into three equation sets: A/B/C, D–J, K/L.
#'     \strong{Use when working with pre-BEC BC inventory data} (Forest
#'     Resources Inventory data from roughly before 2000) or any dataset
#'     where FIZ zone codes are available but BEC zones are not.  Also
#'     ensures reproducibility when replicating older MoF volume
#'     compilations.}
#' }
#'
#' Volume is integrated numerically in 10 cm height slices using Smalian's
#' formula, matching the approach in FAIBBase's \code{treeVolCalculator.R}.
#' Stump volume is computed as a cylinder.
#'
#' @details
#' \strong{KBEC equation} (Kozak 2002/2004):
#' \deqn{d_i = ERR \cdot B_1 \cdot D^{B_2} \cdot H^{B_3} \cdot
#'   X^{(B_4 Z^4 + B_5 e^{-D/H} + B_6 X^{0.1} + B_7/D + B_8 H^Q + B_9 X)}}
#' where \eqn{Z = h_i/H}, \eqn{P = 1.3/H}, \eqn{Q = 1 - Z^{1/3}},
#' \eqn{X = Q/(1 - P^{1/3})}.
#'
#' \strong{KFIZ3 equation} (Kozak, sqrt-transform form):
#' \deqn{d_i = A_0 \cdot D^{A_1} \cdot A_2^D \cdot X^{E_X}}
#' where \eqn{X = (1 - \sqrt{Z})/(1 - \sqrt{P})} (\eqn{P} is a fitted
#' inflection parameter, not fixed at \eqn{1.3/H}),
#' \eqn{E_X = A_3 Z^2 + A_4 \ln(Z+0.001) + A_5\sqrt{X} + A_6(D/H) + A_7 e^Z}.
#'
#' Both equations produce zero DIB at the tree tip and are undefined above
#' total height; \code{btop_height} limits integration for broken-top trees.
#'
#' \strong{Supported species (SP0) for both equations}:
#' \code{"AC"} (poplar/cottonwood), \code{"AT"} (trembling aspen),
#' \code{"B"} (true firs), \code{"C"} (western redcedar),
#' \code{"D"} (red alder), \code{"E"} (paper birch),
#' \code{"F"} (Douglas-fir), \code{"H"} (western hemlock),
#' \code{"L"} (western larch), \code{"MB"} (bigleaf maple),
#' \code{"PA"} (whitebark pine), \code{"PL"} (lodgepole pine),
#' \code{"PW"} (western white pine), \code{"PY"} (ponderosa pine),
#' \code{"S"} (spruce), \code{"Y"} (yellow-cedar).
#'
#' \strong{BEC zones} (for \code{taper_eq = "KBEC"}):
#' \code{"AT"}, \code{"BWBS"}, \code{"CDF"}, \code{"CWH"}, \code{"ESSF"},
#' \code{"ICH"}, \code{"IDF"}, \code{"MH"}, \code{"MS"}, \code{"PP"},
#' \code{"SBPS"}, \code{"SBS"}, \code{"SWB"}.
#'
#' \strong{FIZ zones} (for \code{taper_eq = "KFIZ3"}):
#' \code{"A"} through \code{"L"} (excluding \code{"I"}) -- 12 zones total.
#' Zone groupings: A/B/C = northern/interior, D–J = central/southern,
#' K/L = coastal.
#'
#' @param bec_zone character vector. BC BEC zone code for each tree
#'   (when \code{taper_eq = "KBEC"}), or BC FIZ zone code A–L
#'   (when \code{taper_eq = "KFIZ3"}).  May be scalar (recycled).
#' @param species character vector. BC species code (SP0) for each tree.
#'   May be scalar (recycled to the length of \code{dbh}).
#' @param dbh numeric vector. Diameter at breast height (cm). Must be > 1 cm.
#' @param height numeric vector. Total tree height (m). Must be > 1.4 m.
#' @param volume_type character. One of \code{"WSV"} (whole-stem volume,
#'   default), \code{"MER"} (gross merchantable volume: above stump to the
#'   last height where inside-bark diameter \eqn{\geq} \code{utop_dib}),
#'   \code{"NMR"} (non-merchantable volume: above the merchantable top to the
#'   tree tip or \code{btop_height}; equals \eqn{WSV - STUMP - MER}), or
#'   \code{"STUMP"}.
#' @param stump_height numeric. Stump height (m). Default \code{0.3}.
#'   Matches FAIBBase default.
#' @param utop_dib numeric. Minimum inside-bark diameter (cm) defining the
#'   merchantable top.  Default \code{10} cm.  Only used for \code{"MER"}
#'   and \code{"NMR"}. Matches FAIBBase default.
#' @param btop_height numeric vector or scalar. Height at broken top (m).
#'   Use \code{NA} (default) for sound trees.
#' @param taper_eq character. Taper equation system: \code{"KBEC"} (default)
#'   for modern BEC-stratified data, or \code{"KFIZ3"} for older
#'   FIZ-stratified data.  See Details.
#'
#' @return numeric vector of tree volumes (m\ifelse{html}{\out{<sup>3</sup>}}{$^3$}),
#'   same length as \code{dbh}.  Returns \code{NA} with a warning for trees
#'   whose zone/species combination is not in the coefficient table.
#'
#' @references
#'   Kozak, A. (2004). My last words on taper equations. \emph{The Forestry
#'   Chronicle}, 80(4), 507--515. \doi{10.5558/tfc80507-4}
#'
#'   Luo, Y. (FAIBBase R package, bcgov/FAIBBase).
#'   \code{R/DIB_ICalculator.R} and \code{R/treeVolCalculator.R}.
#'
#' @seealso \code{\link{taper_coefs_kbec}}, \code{\link{taper_coefs_kfiz3}},
#'   \code{\link{volume_citations}}
#' @export
#' @examples
#' # Whole-stem volume, single tree (KBEC, the default)
#' tree_volume(bec_zone = "CWH", species = "H", dbh = 30.7, height = 27.4)
#'
#' # Merchantable volume, multiple trees
#' tree_volume(
#'   bec_zone    = c("CWH", "CWH", "IDF"),
#'   species     = c("H",   "S",   "D"),
#'   dbh         = c(30.7,  42.3,  25.0),
#'   height      = c(27.4,  37.3,  22.0),
#'   volume_type = "MER"
#' )
#'
#' # Same tree using KFIZ3 (FIZ zone required instead of BEC zone)
#' tree_volume(bec_zone = "K", species = "H", dbh = 30.7, height = 27.4,
#'             taper_eq = "KFIZ3")
#'
#' # Scalar zone recycled across trees of the same zone
#' tree_volume("SBS", "PL", dbh = c(15, 20, 25, 30), height = c(12, 16, 20, 24))
tree_volume <- function(bec_zone, species, dbh, height,
                        volume_type  = "WSV",
                        stump_height = 0.3,
                        utop_dib     = 10,
                        btop_height  = NA,
                        taper_eq     = "KBEC") {

  volume_type <- match.arg(volume_type, c("WSV", "MER", "NMR", "STUMP"))
  taper_eq    <- match.arg(taper_eq,    c("KBEC", "KFIZ3"))

  n <- length(dbh)
  if (length(species)     == 1L) species     <- rep(species,     n)
  if (length(bec_zone)    == 1L) bec_zone    <- rep(bec_zone,    n)
  if (length(btop_height) == 1L) btop_height <- rep(btop_height, n)

  if (length(species)  != n) stop("'species' must be scalar or same length as 'dbh'.")
  if (length(bec_zone) != n) stop("'bec_zone' must be scalar or same length as 'dbh'.")
  if (length(height)   != n) stop("'height' must be same length as 'dbh'.")

  result <- numeric(n)

  for (i in seq_len(n)) {
    d  <- dbh[i]; h  <- height[i]
    sp <- species[i]; zone <- bec_zone[i]
    bh <- btop_height[i]

    if (is.na(d) || is.na(h) || d <= 1 || h <= 1.4) {
      warning(sprintf(
        "Tree %d: invalid DBH (%.1f) or height (%.1f); volume set to NA.", i, d, h))
      result[i] <- NA_real_
      next
    }

    if (taper_eq == "KBEC") {
      row <- which(taper_coefs_kbec$species  == sp &
                   taper_coefs_kbec$bec_zone == zone)
      if (length(row) == 0L) {
        warning(sprintf(
          "Tree %d: no KBEC coefficients for species '%s' / BEC zone '%s'; volume set to NA.",
          i, sp, zone))
        result[i] <- NA_real_; next
      }
      cf <- unlist(taper_coefs_kbec[row[1L],
                     c("B1","B2","B3","B4","B5","B6","B7","B8","B9","ERR")])
      dib_fn <- local({ d_ <- d; h_ <- h; cf_ <- cf
                        function(hi) .kozak_kbec_dib(hi, d_, h_, cf_) })
    } else {
      row <- which(taper_coefs_kfiz3$species  == sp &
                   taper_coefs_kfiz3$fiz_zone == zone)
      if (length(row) == 0L) {
        warning(sprintf(
          "Tree %d: no KFIZ3 coefficients for species '%s' / FIZ zone '%s'; volume set to NA.",
          i, sp, zone))
        result[i] <- NA_real_; next
      }
      cf <- unlist(taper_coefs_kfiz3[row[1L],
                     c("P","A0","A1","A2","A3","A4","A5","A6","A7")])
      dib_fn <- local({ d_ <- d; h_ <- h; cf_ <- cf
                        function(hi) .kozak_kfiz3_dib(hi, d_, h_, cf_) })
    }

    result[i] <- .tree_vol_single(d, h, volume_type, stump_height,
                                  utop_dib, bh, dib_fn)
  }
  result
}


# Internal log-section volume calculator -------------------------------------

#' @keywords internal
.tree_vol_sections_single <- function(dbh, height, log_lengths,
                                       stump_height, utop_dib,
                                       btop_height, dib_fn) {
  vcons <- pi / 40000

  # Effective upper limit (broken top or tree tip)
  upper <- if (!is.na(btop_height) && btop_height < height) btop_height else height

  # Section boundaries above the stump
  bots <- stump_height + c(0, cumsum(log_lengths[-length(log_lengths)]))
  tops <- stump_height + cumsum(log_lengths)

  # Keep only sections that start below the effective upper limit
  in_range <- bots < upper
  if (!any(in_range)) {
    return(data.frame(section = integer(0), height_bot = numeric(0),
                      height_top = numeric(0), length = numeric(0),
                      dib_bot = numeric(0), dib_top = numeric(0),
                      vol_gross = numeric(0), vol_merch = numeric(0),
                      stringsAsFactors = FALSE))
  }
  bots <- bots[in_range]
  tops <- pmin(tops[in_range], upper)
  n    <- length(bots)

  # Fine-grained height grid (10 cm steps) with section boundaries forced in
  grid <- seq(stump_height, upper, by = 0.1)
  if (tail(grid, 1) < upper) grid <- c(grid, upper)
  grid <- sort(unique(round(c(grid, bots, tops), 10)))
  grid <- grid[grid >= stump_height & grid <= upper]

  dibs <- dib_fn(grid)
  dib_bh_v <- dib_fn(1.3)
  below_bh  <- grid < 1.3
  if (any(below_bh)) dibs[below_bh] <- pmax(dibs[below_bh], dib_bh_v)
  dibs[!is.finite(dibs) | dibs < 0] <- 0

  areas      <- vcons * dibs^2
  lens       <- diff(grid)
  vol_slices <- (head(areas, -1) + tail(areas, -1)) / 2 * lens
  slice_tops     <- grid[-1]
  slice_dib_tops <- dibs[-1]  # DIB at top of each slice

  # Assign each slice to a section: section k spans (bots[k], tops[k]]
  breaks      <- c(bots[1] - 1e-9, tops)
  section_idx <- as.integer(cut(slice_tops, breaks = breaks, labels = FALSE))

  results <- vector("list", n)
  for (k in seq_len(n)) {
    in_sec  <- !is.na(section_idx) & section_idx == k
    v_gross <- sum(vol_slices[in_sec])
    v_merch <- sum(vol_slices[in_sec & slice_dib_tops >= utop_dib])

    dib_b <- dib_fn(bots[k])
    dib_t <- dib_fn(tops[k])
    dib_b <- max(0, if (is.finite(dib_b)) dib_b else 0)
    dib_t <- max(0, if (is.finite(dib_t)) dib_t else 0)

    results[[k]] <- data.frame(
      section    = k,
      height_bot = bots[k],
      height_top = tops[k],
      length     = tops[k] - bots[k],
      dib_bot    = round(dib_b, 3),
      dib_top    = round(dib_t, 3),
      vol_gross  = v_gross,
      vol_merch  = v_merch,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, results)
}


# Log-section volume (exported) -----------------------------------------------

#' Calculate log-section volumes for a single tree
#'
#' @description
#' Partitions a tree stem into log sections and returns the gross and
#' merchantable inside-bark volume
#' (m\ifelse{html}{\out{<sup>3</sup>}}{$^3$}) for each section.
#' Supports both KBEC and KFIZ3 taper equations; see
#' \code{\link{tree_volume}} Details for guidance on which to use.
#'
#' Sections are stacked above the stump: section 1 runs from
#' \code{stump_height} to \code{stump_height + log_lengths[1]}, section 2
#' continues from there, and so on.  Any section extending above the
#' effective tree top (or \code{btop_height}) is truncated.
#'
#' This function operates on a \strong{single tree}.  For multiple trees
#' use \code{Map()} or \code{lapply()} -- see Examples.
#'
#' @inheritParams tree_volume
#' @param log_lengths numeric vector. Length (m) of each log section in
#'   ascending order from stump to top. All values must be positive.
#'
#' @return A data frame with one row per section:
#' \describe{
#'   \item{\code{section}}{Integer section index (1 = lowest).}
#'   \item{\code{height_bot}}{Bottom height of the section (m).}
#'   \item{\code{height_top}}{Top height of the section (m).}
#'   \item{\code{length}}{Effective section length after truncation (m).}
#'   \item{\code{dib_bot}}{Inside-bark diameter at section bottom (cm).}
#'   \item{\code{dib_top}}{Inside-bark diameter at section top (cm).}
#'   \item{\code{vol_gross}}{Gross inside-bark volume
#'     (m\ifelse{html}{\out{<sup>3</sup>}}{$^3$}); sum of 10 cm Smalian
#'     slices within the section.}
#'   \item{\code{vol_merch}}{Merchantable inside-bark volume
#'     (m\ifelse{html}{\out{<sup>3</sup>}}{$^3$}); slices where
#'     slice-top DIB \eqn{\geq} \code{utop_dib}.}
#' }
#' Returns \code{NULL} with a warning for invalid inputs or unknown
#' zone/species combinations.
#'
#' @export
#' @seealso \code{\link{tree_volume}}, \code{\link{tree_profile}},
#'   \code{\link{volume_citations}}
#' @references
#'   Kozak, A. (2004). My last words on taper equations. \emph{The Forestry
#'   Chronicle}, 80(4), 507--515. \doi{10.5558/tfc80507-4}
#'
#'   Luo, Y. (FAIBBase R package, bcgov/FAIBBase).
#'   \code{R/treeVolCalculator.R} -- log-section volume implementation.
#' @examples
#' # Four 3-m sections, KBEC (default)
#' tree_volume_section("CWH", "H", dbh = 30.7, height = 27.4,
#'                     log_lengths = rep(3, 4))
#'
#' # Same tree, KFIZ3 (coastal FIZ zone K)
#' tree_volume_section("K", "H", dbh = 30.7, height = 27.4,
#'                     log_lengths = rep(3, 4), taper_eq = "KFIZ3")
#'
#' # Section gross volumes should sum to approximately WSV minus stump volume
#' secs <- tree_volume_section("CWH", "H", 30.7, 27.4, log_lengths = rep(3, 9))
#' cat("Sum section vol_gross :", round(sum(secs$vol_gross), 5), "m3\n")
#' cat("WSV minus STUMP       :",
#'     round(tree_volume("CWH", "H", 30.7, 27.4, "WSV") -
#'           tree_volume("CWH", "H", 30.7, 27.4, "STUMP"), 5), "m3\n")
#'
#' # Multiple trees using Map()
#' \dontrun{
#' dat  <- data.frame(BEC = c("CWH", "SBS"), SP = c("H", "PL"),
#'                    DBH = c(30.7, 22.0), HT = c(27.4, 18.0))
#' secs <- Map(tree_volume_section,
#'             bec_zone = dat$BEC, species = dat$SP,
#'             dbh = dat$DBH, height = dat$HT,
#'             log_lengths = list(rep(3, 4), rep(2.5, 4)))
#' }
tree_volume_section <- function(bec_zone, species, dbh, height,
                                log_lengths,
                                stump_height = 0.3,
                                utop_dib     = 10,
                                btop_height  = NA,
                                taper_eq     = "KBEC") {
  if (length(dbh) != 1L || length(height) != 1L)
    stop("'tree_volume_section' operates on a single tree; use Map() for multiple trees.")
  if (!is.numeric(log_lengths) || length(log_lengths) == 0L || any(log_lengths <= 0))
    stop("'log_lengths' must be a non-empty numeric vector of positive values.")
  taper_eq <- match.arg(taper_eq, c("KBEC", "KFIZ3"))

  if (is.na(dbh) || !is.finite(dbh) || dbh <= 1) {
    warning("tree_volume_section: invalid DBH (<= 1 cm). Returning NULL.")
    return(NULL)
  }
  if (is.na(height) || !is.finite(height) || height <= 1.4) {
    warning("tree_volume_section: invalid height (<= 1.4 m). Returning NULL.")
    return(NULL)
  }

  if (taper_eq == "KBEC") {
    row <- which(taper_coefs_kbec$species  == species &
                 taper_coefs_kbec$bec_zone == bec_zone)
    if (length(row) == 0L) {
      warning(sprintf(
        "tree_volume_section: no KBEC coefficients for species '%s' / BEC zone '%s'. Returning NULL.",
        species, bec_zone))
      return(NULL)
    }
    cf <- unlist(taper_coefs_kbec[row[1L],
                   c("B1","B2","B3","B4","B5","B6","B7","B8","B9","ERR")])
    dib_fn <- local({ d_ <- dbh; h_ <- height; cf_ <- cf
                      function(hi) .kozak_kbec_dib(hi, d_, h_, cf_) })
  } else {
    row <- which(taper_coefs_kfiz3$species  == species &
                 taper_coefs_kfiz3$fiz_zone == bec_zone)
    if (length(row) == 0L) {
      warning(sprintf(
        "tree_volume_section: no KFIZ3 coefficients for species '%s' / FIZ zone '%s'. Returning NULL.",
        species, bec_zone))
      return(NULL)
    }
    cf <- unlist(taper_coefs_kfiz3[row[1L],
                   c("P","A0","A1","A2","A3","A4","A5","A6","A7")])
    dib_fn <- local({ d_ <- dbh; h_ <- height; cf_ <- cf
                      function(hi) .kozak_kfiz3_dib(hi, d_, h_, cf_) })
  }

  .tree_vol_sections_single(dbh, height, log_lengths,
                             stump_height, utop_dib, btop_height, dib_fn)
}


# Internal stem-profile builder -----------------------------------------------

#' @keywords internal
.tree_profile_build <- function(dbh, height, stump_height, utop_dib,
                                btop_height, dib_fn) {
  vcons <- pi / 40000  # DIB (cm)^2 x height (m) -> m^3

  # Effective stump DIB (floor at breast-height DIB, matching the compiler)
  dib_stump <- dib_fn(stump_height)
  dib_bh    <- dib_fn(1.3)
  if (dib_bh > dib_stump) dib_stump <- dib_bh

  # Upper boundary: break point for broken-top trees, tree tip otherwise
  upper <- if (!is.na(btop_height) && btop_height < height) btop_height else height

  # Height grid at 10 cm steps; force 1.3 m in for the breast-height annotation
  hts <- round(seq(stump_height, upper, by = 0.1), 1)
  if (1.3 > stump_height && 1.3 < upper) {
    hts <- sort(unique(c(hts, 1.3)))
  }
  if (tail(hts, 1) < upper) hts <- c(hts, upper)

  dibs <- dib_fn(hts)
  # Apply breast-height floor (matches FAIBBase treeProfile behaviour)
  below_bh_p <- hts < 1.3
  if (any(below_bh_p)) dibs[below_bh_p] <- pmax(dibs[below_bh_p], dib_bh)
  dibs[!is.finite(dibs) | dibs < 0] <- 0

  # Smalian slice volumes; last row has no "slice above" so vol_slice = 0
  areas  <- vcons * dibs^2
  lens   <- diff(hts)
  vslice <- c((head(areas, -1) + tail(areas, -1)) / 2 * lens, 0)

  # Prepend the stump-cylinder row at height = 0
  vstump <- vcons * stump_height * dib_stump^2
  hts    <- c(0, hts)
  dibs   <- c(dib_stump, dibs)
  vslice <- c(vstump, vslice)
  cumvol <- cumsum(vslice)

  # Annotations
  comment <- rep("", length(hts))
  tol <- 1e-6
  comment[abs(hts - 0)            < tol] <- "ground"
  comment[abs(hts - stump_height) < tol] <- "stump height"
  # breast height (may coincide with stump height for very short trees)
  bh_idx <- which(abs(hts - 1.3) < tol)
  if (length(bh_idx) > 0L) comment[bh_idx] <- "breast height"

  # Max merchantable height: last row where dib >= utop_dib (above stump)
  mer_mask <- dibs >= utop_dib & hts > stump_height
  if (any(mer_mask)) {
    mer_idx <- max(which(mer_mask))
    comment[mer_idx] <- trimws(paste(comment[mer_idx], "max merchantable height"))
  }

  # Tip or break-point label
  n <- length(hts)
  tip_label <- if (!is.na(btop_height) && btop_height < height) "break height" else "tip"
  comment[n] <- trimws(paste(comment[n], tip_label))

  data.frame(height = hts, dib = dibs, vol_slice = vslice,
             cumvol = cumvol, comment = comment,
             stringsAsFactors = FALSE)
}


# Exported stem-profile function ----------------------------------------------

#' Generate a stem taper profile for a single tree
#'
#' @description
#' Returns a data frame of inside-bark diameter (DIB) and cumulative volume at
#' every 10 cm height step along the stem.  Supports both KBEC and KFIZ3
#' taper equations; see \code{\link{tree_volume}} Details for guidance on
#' which to use.
#'
#' The profile is useful for:
#' \itemize{
#'   \item \strong{Log bucking / grade simulation} -- knowing DIB at each
#'     10 cm lets you locate where any minimum-diameter threshold is met.
#'   \item \strong{Custom volume segments} -- subtract \code{cumvol} at two
#'     heights to get the volume of any bolt or section.
#'   \item \strong{Taper model validation} -- compare predicted DIB against
#'     independently measured upper-stem diameters.
#' }
#'
#' @details
#' The profile begins at height = 0 (the stump-cylinder base row), followed by
#' rows at 0.1 m intervals from \code{stump_height} to total height (or
#' \code{btop_height} for broken-top trees).  Breast height (1.3 m) is always
#' included in the grid.
#'
#' The \code{comment} column marks key heights:
#' \describe{
#'   \item{\code{"ground"}}{height = 0; the stump-cylinder base}
#'   \item{\code{"stump height"}}{height = \code{stump_height}}
#'   \item{\code{"breast height"}}{height = 1.3 m}
#'   \item{\code{"max merchantable height"}}{last height where DIB >=
#'     \code{utop_dib}}
#'   \item{\code{"break height"}}{the stem break point for broken-top trees}
#'   \item{\code{"tip"}}{the tree tip (non-broken-top trees)}
#' }
#'
#' This function operates on a **single tree**.  To generate profiles for many
#' trees use \code{Map()} or \code{lapply()} — see Examples.
#'
#' @inheritParams tree_volume
#'
#' @return A data frame with five columns:
#' \describe{
#'   \item{\code{height}}{Height above ground (m).}
#'   \item{\code{dib}}{Inside-bark diameter at this height (cm).}
#'   \item{\code{vol_slice}}{Volume (m³) from this height to the next,
#'     computed by Smalian's formula.  Zero for the final row.}
#'   \item{\code{cumvol}}{Cumulative volume (m³) from ground to this height.}
#'   \item{\code{comment}}{Annotation; see Details.}
#' }
#' Returns \code{NULL} with a warning for invalid inputs or unknown
#' species/BEC combinations.
#'
#' @inheritParams tree_volume
#' @export
#' @seealso \code{\link{tree_volume}}, \code{\link{taper_coefs_kbec}},
#'   \code{\link{taper_coefs_kfiz3}}
#' @references
#'   Kozak, A. (2004). My last words on taper equations.
#'   \emph{The Forestry Chronicle}, 80(4), 507--515.
#'   \doi{10.5558/tfc80507-4}
#'
#' @examples
#' # Full stem profile (KBEC, default)
#' prof <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
#' head(prof)
#'
#' # Same tree using KFIZ3 (coastal FIZ zone K)
#' prof_fiz <- tree_profile("K", "H", dbh = 30.7, height = 27.4,
#'                          taper_eq = "KFIZ3")
#'
#' # Volume from ground to 10 m
#' prof$cumvol[which.min(abs(prof$height - 10))]
#'
#' # Multiple trees: use Map()
#' \dontrun{
#' profiles <- Map(tree_profile,
#'   bec_zone = df$BEC_ZONE, species = df$SPECIES,
#'   dbh = df$DBH, height = df$HEIGHT)
#' }
tree_profile <- function(bec_zone, species, dbh, height,
                         stump_height = 0.3,
                         utop_dib     = 10,
                         btop_height  = NA,
                         taper_eq     = "KBEC") {
  if (length(dbh) != 1L || length(height) != 1L)
    stop("'tree_profile' operates on a single tree; 'dbh' and 'height' must be scalar.")
  taper_eq <- match.arg(taper_eq, c("KBEC", "KFIZ3"))

  if (!is.finite(dbh) || dbh <= 1) {
    warning("tree_profile: invalid DBH (<= 1 cm). Returning NULL.")
    return(NULL)
  }
  if (!is.finite(height) || height <= 1.4) {
    warning("tree_profile: invalid height (<= 1.4 m). Returning NULL.")
    return(NULL)
  }

  if (taper_eq == "KBEC") {
    row <- which(taper_coefs_kbec$species  == species &
                 taper_coefs_kbec$bec_zone == bec_zone)
    if (length(row) == 0L) {
      warning(sprintf(
        "tree_profile: no KBEC coefficients for species '%s' / BEC zone '%s'. Returning NULL.",
        species, bec_zone))
      return(NULL)
    }
    cf <- unlist(taper_coefs_kbec[row[1L],
                   c("B1","B2","B3","B4","B5","B6","B7","B8","B9","ERR")])
    dib_fn <- local({ d_ <- dbh; h_ <- height; cf_ <- cf
                      function(hi) .kozak_kbec_dib(hi, d_, h_, cf_) })
  } else {
    row <- which(taper_coefs_kfiz3$species  == species &
                 taper_coefs_kfiz3$fiz_zone == bec_zone)
    if (length(row) == 0L) {
      warning(sprintf(
        "tree_profile: no KFIZ3 coefficients for species '%s' / FIZ zone '%s'. Returning NULL.",
        species, bec_zone))
      return(NULL)
    }
    cf <- unlist(taper_coefs_kfiz3[row[1L],
                   c("P","A0","A1","A2","A3","A4","A5","A6","A7")])
    dib_fn <- local({ d_ <- dbh; h_ <- height; cf_ <- cf
                      function(hi) .kozak_kfiz3_dib(hi, d_, h_, cf_) })
  }

  .tree_profile_build(dbh, height, stump_height, utop_dib, btop_height, dib_fn)
}


# Back-calculate total height from a broken-top measurement ------------------

#' Back-calculate total tree height from a broken-top measurement
#'
#' @description
#' Estimates the total height a broken-top tree would have reached had it
#' remained intact, using a Kozak taper equation.  Total height is found by
#' bisection (via \code{\link[stats]{uniroot}}): the search identifies the
#' total height \eqn{H} at which the taper equation predicts an inside-bark
#' diameter equal to \code{btop_dib} at \code{btop_height}.
#'
#' This mirrors the approach of \code{FAIBBase::heightEstimateForBTOP_D()},
#' which is used when field crews recorded the break height and the inside-bark
#' diameter at the break rather than a visually projected total height.
#'
#' @details
#' The bisection search spans \code{(btop_height, 66.6]} metres; 66.6 m is
#' the maximum plausible tree height used in FAIBBase.  Returns \code{NA} with
#' a warning when:
#' \itemize{
#'   \item any required input is \code{NA};
#'   \item \code{btop_dib >= dbh} (DIB at break cannot exceed DBH);
#'   \item \code{btop_height < 1.4} or \code{btop_height > 60} (outside
#'     valid range, matching FAIBBase guards);
#'   \item no taper coefficients exist for the species/zone combination; or
#'   \item no root is found in the search interval (unusual: would require an
#'     implausibly large tree).
#' }
#'
#' @param dbh numeric vector. Diameter at breast height (cm).
#' @param btop_height numeric vector. Height of the broken top above ground
#'   (m), i.e., the stump height of the broken section.
#' @param btop_dib numeric vector. Inside-bark diameter (cm) measured at
#'   \code{btop_height}.
#' @param bec_zone character vector. BC BEC zone code (when
#'   \code{taper_eq = "KBEC"}), or FIZ zone code A--L (when
#'   \code{taper_eq = "KFIZ3"}).  May be scalar (recycled).
#' @param species character vector. BC species code (SP0).  May be scalar.
#' @param taper_eq character. Taper equation system: \code{"KBEC"} (default)
#'   or \code{"KFIZ3"}.  See \code{\link{tree_volume}} for guidance.
#'
#' @return numeric vector of estimated total heights (m), rounded to 0.1 m.
#'   Same length as \code{dbh}.  \code{NA} is returned for trees that fail
#'   validation or where no solution exists.
#'
#' @references
#'   Kozak, A. (2004). My last words on taper equations. \emph{The Forestry
#'   Chronicle}, 80(4), 507--515. \doi{10.5558/tfc80507-4}
#'
#'   Luo, Y. (FAIBBase R package, bcgov/FAIBBase).
#'   \code{R/heightEstimateForBTOP.R} (\code{heightEstimateForBTOP_D}).
#'
#' @seealso \code{\link{tree_volume}}, \code{\link{tree_profile}},
#'   \code{\link{volume_citations}}
#' @importFrom stats uniroot
#' @export
#' @examples
#' # Single broken-top tree in the CWH zone
#' ht_from_btop(dbh = 42.0, btop_height = 22.5, btop_dib = 8.3,
#'              bec_zone = "CWH", species = "H")
#'
#' # Multiple trees
#' ht_from_btop(
#'   dbh        = c(42.0, 30.5, 55.2),
#'   btop_height = c(22.5, 18.0, 30.1),
#'   btop_dib   = c( 8.3,  6.1, 10.4),
#'   bec_zone   = "SBS",
#'   species    = c("S", "PL", "S")
#' )
ht_from_btop <- function(dbh, btop_height, btop_dib,
                          bec_zone, species,
                          taper_eq = "KBEC") {

  taper_eq <- match.arg(taper_eq, c("KBEC", "KFIZ3"))

  n <- length(dbh)
  if (length(species)     == 1L) species     <- rep(species,     n)
  if (length(bec_zone)    == 1L) bec_zone    <- rep(bec_zone,    n)
  if (length(btop_height) == 1L) btop_height <- rep(btop_height, n)
  if (length(btop_dib)    == 1L) btop_dib    <- rep(btop_dib,    n)

  result <- numeric(n)

  for (i in seq_len(n)) {
    d    <- dbh[i]
    hb   <- btop_height[i]
    dib  <- btop_dib[i]
    sp   <- species[i]
    zone <- bec_zone[i]

    # --- Input validation (mirrors FAIBBase guards) -------------------------
    if (anyNA(c(d, hb, dib, zone, sp))) {
      result[i] <- NA_real_; next
    }
    if (dib >= d) {
      warning(sprintf(
        "Tree %d: btop_dib (%.2f) >= dbh (%.2f); cannot solve. Returning NA.", i, dib, d))
      result[i] <- NA_real_; next
    }
    if (hb < 1.4 || hb > 60) {
      warning(sprintf(
        "Tree %d: btop_height (%.1f) outside valid range [1.4, 60] m. Returning NA.", i, hb))
      result[i] <- NA_real_; next
    }

    # --- Coefficient lookup ------------------------------------------------
    if (taper_eq == "KBEC") {
      row <- which(taper_coefs_kbec$species  == sp &
                   taper_coefs_kbec$bec_zone == zone)
      if (length(row) == 0L) {
        warning(sprintf(
          "Tree %d: no KBEC coefficients for species '%s' / BEC zone '%s'. Returning NA.",
          i, sp, zone))
        result[i] <- NA_real_; next
      }
      cf <- unlist(taper_coefs_kbec[row[1L],
                     c("B1","B2","B3","B4","B5","B6","B7","B8","B9","ERR")])
      dib_at <- function(H) .kozak_kbec_dib(hb, d, H, cf)
    } else {
      row <- which(taper_coefs_kfiz3$species  == sp &
                   taper_coefs_kfiz3$fiz_zone == zone)
      if (length(row) == 0L) {
        warning(sprintf(
          "Tree %d: no KFIZ3 coefficients for species '%s' / FIZ zone '%s'. Returning NA.",
          i, sp, zone))
        result[i] <- NA_real_; next
      }
      cf <- unlist(taper_coefs_kfiz3[row[1L],
                     c("P","A0","A1","A2","A3","A4","A5","A6","A7")])
      dib_at <- function(H) .kozak_kfiz3_dib(hb, d, H, cf)
    }

    # --- Bisection via uniroot ---------------------------------------------
    # Objective: find H in (hb, 66.6] where taper-predicted DIB at hb == dib
    H_min <- hb + 0.01
    H_max <- 66.6

    f_min <- dib_at(H_min) - dib
    f_max <- dib_at(H_max) - dib

    if (is.na(f_min) || is.na(f_max) || f_min * f_max > 0) {
      warning(sprintf(
        "Tree %d: no solution in [%.1f, 66.6] m (btop_dib may be implausible). Returning NA.",
        i, H_min))
      result[i] <- NA_real_; next
    }

    sol <- tryCatch(
      uniroot(function(H) dib_at(H) - dib,
              lower = H_min, upper = H_max, tol = 0.01)$root,
      error = function(e) NA_real_
    )
    result[i] <- if (is.na(sol)) NA_real_ else round(sol, 1L)
  }
  result
}
