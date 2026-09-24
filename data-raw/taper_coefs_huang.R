# data-raw/taper_coefs_huang.R
#
# Builds the taper_coefs_huang package dataset.
#
# Source:
#   Huang, S. (1994). Ecologically based individual tree volume estimation
#   for major Alberta tree species. Alberta Environmental Protection, Land
#   and Forest Services, Forest Management Division. (Report 1: Individual
#   tree volume estimation procedures for Alberta.)
#
#   Coefficients transcribed from the CTAE R package (ptompalski/CTAE),
#   object parameters_HuangV, which digitises Huang et al. (1994) Appendix 3.
#
# The Huang variable-exponent taper equation (Kozak sqrt-transform form):
#   dib = a0 * DBH^a1 * a2^DBH *
#         X^(b1*Z^2 + b2*log(Z+0.001) + b3*sqrt(Z) + b4*exp(Z) + b5*DBH/H)
# where:
#   Z = h_i / H                          (relative height)
#   X = (1 - sqrt(Z)) / (1 - sqrt(0.225))
#
# The long CSV (one row per species x subregion x parameter) is pivoted to a
# wide table: one row per species x natural subregion, with columns
# a0, a1, a2, b1, b2, b3, b4, b5.

library(usethis)

raw <- read.csv("data-raw/parameters_HuangV.csv", stringsAsFactors = FALSE)

# Drop the leading unnamed row-number column if present
raw <- raw[, c("species", "parameter", "estimate",
               "NaturalSubregionNum", "NaturalSubregionCode")]

# Pivot parameters (a0..b5) from long to wide
taper_coefs_huang <- reshape(
  raw,
  idvar     = c("species", "NaturalSubregionNum", "NaturalSubregionCode"),
  timevar   = "parameter",
  direction = "wide"
)
names(taper_coefs_huang) <- sub("^estimate\\.", "", names(taper_coefs_huang))

# Order columns and rows
taper_coefs_huang <- taper_coefs_huang[, c(
  "species", "NaturalSubregionCode", "NaturalSubregionNum",
  "a0", "a1", "a2", "b1", "b2", "b3", "b4", "b5"
)]
taper_coefs_huang <- taper_coefs_huang[
  order(taper_coefs_huang$species, taper_coefs_huang$NaturalSubregionCode), ]
rownames(taper_coefs_huang) <- NULL

# Only keep rows with a complete set of 8 coefficients
complete <- stats::complete.cases(
  taper_coefs_huang[, c("a0", "a1", "a2", "b1", "b2", "b3", "b4", "b5")])
taper_coefs_huang <- taper_coefs_huang[complete, ]
rownames(taper_coefs_huang) <- NULL

usethis::use_data(taper_coefs_huang, overwrite = TRUE)
