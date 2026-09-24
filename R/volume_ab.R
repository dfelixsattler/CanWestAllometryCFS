# Huang (1994) variable-exponent taper / volume equations for Alberta tree
# species.  Provides tree_volume_huang() for individual-tree merchantable and
# whole-stem volume, stratified by species and Alberta natural subregion.
#
# Exported function: tree_volume_huang().

# Internal Huang inside-bark diameter -----------------------------------------

#' @keywords internal
.huang_dib <- function(height_i, dbh, height, cf) {
  # Inside-bark diameter (cm) at height_i (m) for a single tree.
  # cf: named numeric vector with a0, a1, a2, b1, b2, b3, b4, b5.
  z  <- height_i / height
  x  <- (1 - sqrt(z)) / (1 - sqrt(0.225))
  ex <- cf[["b1"]] * z^2 +
        cf[["b2"]] * log(z + 0.001) +
        cf[["b3"]] * sqrt(z) +
        cf[["b4"]] * exp(z) +
        cf[["b5"]] * dbh / height
  dib <- (cf[["a0"]] * dbh^cf[["a1"]]) * (cf[["a2"]]^dbh) * x^ex
  dib[!is.finite(dib) | x <= 0] <- 0
  dib[dib < 0] <- 0
  dib
}

# Internal per-tree Huang volume ----------------------------------------------

#' @keywords internal
.huang_vol_single <- function(dbh, height, cf, utop_dib, stump_height) {
  # Returns c(merch, total) for one tree following Huang et al. (1994) Appendix 3.
  if (!is.finite(dbh) || !is.finite(height) || dbh <= 0 || height <= 0) {
    return(c(merch = NA_real_, total = NA_real_))
  }

  # Iteratively solve for g0 = merchantable-height ratio (h_merch / H) at the
  # point where inside-bark diameter equals utop_dib.
  g0 <- 0.9
  g1 <- 0
  iter <- 0L
  repeat {
    if (!is.finite(g0) || !is.finite(g1) ||
        abs(g0 - g1) <= 1e-08 || iter > 1000L) break
    cc <- cf[["b1"]] * g0^2 + cf[["b2"]] * log(g0 + 0.001) +
          cf[["b3"]] * sqrt(g0) + cf[["b4"]] * exp(g0) +
          cf[["b5"]] * (dbh / height)
    base <- utop_dib / (cf[["a0"]] * dbh^cf[["a1"]] * cf[["a2"]]^dbh)
    g1 <- (1 - (base^(1 / cc)) * (1 - sqrt(0.225)))^2
    g0 <- (g0 + g1) / 2
    iter <- iter + 1L
  }

  if (iter > 1000L || !is.finite(g1) || !is.finite(g0)) {
    return(c(merch = NA_real_, total = NA_real_))
  }

  # Merchantable height (hi) and 20 equal sections from the stump to hi
  hi   <- g0 * height
  mlen <- 1:20 * (hi - stump_height) / 20 + stump_height

  # Inside-bark diameter at the stump and at each section boundary
  dib_stump <- .huang_dib(stump_height, dbh, height, cf)
  dib_secs  <- .huang_dib(mlen, dbh, height, cf)
  dibm      <- c(dib_stump, dib_secs)   # 21 values

  # Merchantable volume by Newton's formula over 10 paired sections
  k <- 0.00007854 * (((hi - stump_height) / 10) / 6)
  idx <- seq(1L, 19L, by = 2L)
  mvol <- sum(vapply(idx, function(i)
    k * (dibm[i]^2 + 4 * dibm[i + 1L]^2 + dibm[i + 2L]^2), numeric(1)))

  # Tip cone above merchantable height + stump cylinder
  tipvol <- 0.00007854 * dibm[21]^2 * (height - hi) / 3
  volstp <- 0.00007854 * dibm[1]^2 * stump_height
  tvol   <- mvol + tipvol + volstp

  c(merch = mvol, total = tvol)
}


#' Calculate individual tree volume using the Huang (1994) Alberta taper equation
#'
#' @description
#' Computes individual-tree merchantable and whole-stem volume for major
#' western Canadian tree species using the Huang et al. (1994) ecologically
#' based variable-exponent taper equation, with coefficients stratified by
#' species and Alberta natural subregion.
#'
#' This is a stand-alone western-Canada feature and is independent of the
#' BC \code{\link{tree_volume}} (Kozak KBEC/KFIZ3) equations.  Use this
#' function for Alberta, Saskatchewan, and Manitoba plot data.
#'
#' @details
#' The Huang taper equation (Kozak sqrt-transform form) predicts inside-bark
#' diameter \eqn{d_i} at height \eqn{h_i}:
#' \deqn{d_i = a_0 \cdot D^{a_1} \cdot a_2^{D} \cdot
#'   X^{(b_1 Z^2 + b_2 \ln(Z + 0.001) + b_3 \sqrt{Z} + b_4 e^{Z} + b_5 D/H)}}
#' where \eqn{Z = h_i / H} and \eqn{X = (1 - \sqrt{Z})/(1 - \sqrt{0.225})}.
#'
#' Merchantable volume is integrated by Newton's formula over 20 equal-length
#' sections from the stump height to the merchantable height (the point at
#' which inside-bark diameter falls to \code{utop_dib}).  Whole-stem volume
#' additionally includes the tip cone above merchantable height and the stump
#' cylinder.
#'
#' \strong{Utilization standards.}  The defaults follow a 13 / 7 / 30 standard:
#' a minimum 13 cm DBH, a 7 cm minimum top inside-bark diameter, and a 30 cm
#' (0.3 m) stump.  Trees with \code{dbh < min_dbh} return \code{NA} for both
#' merchantable and whole-stem volume, because they fall below the utilization
#' standard.  Set \code{min_dbh = 0} to obtain volumes for all trees regardless
#' of DBH.
#'
#' @param species  character vector. Huang species code (genus.species form),
#'   e.g. \code{"PINU.CON"} (lodgepole pine), \code{"POPU.TRE"} (trembling
#'   aspen), \code{"PICE.GLA"} (white spruce).  See
#'   \code{unique(taper_coefs_huang$species)} for the full list.  May be scalar
#'   (recycled).
#' @param dbh      numeric vector. Diameter at breast height (cm).
#' @param height   numeric vector. Total tree height (m).
#' @param subregion character vector. Alberta natural subregion code
#'   (e.g. \code{"LF"} Lower Foothills, \code{"CM"} Central Mixedwood).
#'   The Huang model is stratified by subregion, so this must be supplied.
#'   Only black spruce (\code{"PICE.MAR"}) has a \code{"Province"}-wide
#'   fit.  See \code{unique(taper_coefs_huang$NaturalSubregionCode)} for the
#'   available codes.  May be scalar (recycled).
#' @param min_dbh  numeric. Minimum DBH (cm) for a tree to be considered
#'   merchantable.  Trees below this return \code{NA}.  Default \code{13}.
#'   Set to \code{0} to compute volume for all trees.
#' @param utop_dib numeric. Minimum top inside-bark diameter (cm) defining the
#'   merchantable top.  Default \code{7}.
#' @param stump_height numeric. Stump height (m).  Default \code{0.3}.
#'
#' @return A data frame with \code{length(dbh)} rows and two columns:
#'   \describe{
#'     \item{\code{merch}}{Gross merchantable volume
#'       (m\ifelse{html}{\out{<sup>3</sup>}}{\eqn{^3}}) between the stump and
#'       the merchantable top (no tip, no stump).}
#'     \item{\code{total}}{Whole-stem volume
#'       (m\ifelse{html}{\out{<sup>3</sup>}}{\eqn{^3}}) including the tip cone
#'       and stump; equivalent to WSV in the BC compiler.}
#'   }
#'   \code{NA} is returned for trees below \code{min_dbh}, with invalid inputs,
#'   or where the species/subregion combination has no coefficients.
#'
#' @references
#'   Huang, S. (1994). Ecologically based individual tree volume estimation
#'   for major Alberta tree species. Alberta Environmental Protection, Land
#'   and Forest Services, Forest Management Division.
#'
#' @seealso \code{\link{tree_volume}} for BC Kozak taper volume,
#'   \code{\link{taper_coefs_huang}}
#' @export
#' @examples
#' # Single lodgepole pine in the Lower Foothills subregion
#' tree_volume_huang(species = "PINU.CON", dbh = 25, height = 20,
#'                   subregion = "LF")
#'
#' # Multiple trees in the Central Mixedwood subregion
#' tree_volume_huang(
#'   species   = c("PINU.CON", "POPU.TRE", "PICE.GLA"),
#'   dbh       = c(25, 30, 40),
#'   height    = c(20, 24, 28),
#'   subregion = "CM"
#' )
#'
#' # Custom utilization standard: 15 cm DBH, 10 cm top, 15 cm stump
#' tree_volume_huang("PINU.CON", dbh = 25, height = 20, subregion = "LF",
#'                   min_dbh = 15, utop_dib = 10, stump_height = 0.15)
tree_volume_huang <- function(species, dbh, height,
                              subregion,
                              min_dbh      = 13,
                              utop_dib     = 7,
                              stump_height = 0.3) {

  n <- length(dbh)
  if (length(species)   == 1L) species   <- rep(species,   n)
  if (length(subregion) == 1L) subregion <- rep(subregion, n)

  if (length(species)   != n) stop("'species' must be scalar or same length as 'dbh'.")
  if (length(subregion) != n) stop("'subregion' must be scalar or same length as 'dbh'.")
  if (length(height)    != n) stop("'height' must be same length as 'dbh'.")

  # Upper Boreal Highlands (UB) share coefficients with Lower Boreal Highlands
  subregion[subregion == "UB"] <- "LBH"

  merch <- numeric(n)
  total <- numeric(n)

  cf_cols <- c("a0", "a1", "a2", "b1", "b2", "b3", "b4", "b5")

  for (i in seq_len(n)) {
    d  <- dbh[i]; h <- height[i]
    sp <- species[i]; sr <- subregion[i]

    if (is.na(d) || is.na(h) || d < min_dbh) {
      merch[i] <- NA_real_; total[i] <- NA_real_; next
    }

    row <- which(taper_coefs_huang$species              == sp &
                 taper_coefs_huang$NaturalSubregionCode == sr)
    if (length(row) == 0L) {
      warning(sprintf(
        "Tree %d: no Huang coefficients for species '%s' / subregion '%s'; volume set to NA.",
        i, sp, sr))
      merch[i] <- NA_real_; total[i] <- NA_real_; next
    }

    cf <- unlist(taper_coefs_huang[row[1L], cf_cols])
    v  <- .huang_vol_single(d, h, cf, utop_dib, stump_height)
    merch[i] <- v[["merch"]]
    total[i] <- v[["total"]]
  }

  data.frame(merch = merch, total = total)
}
