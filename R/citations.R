# Citation helpers that return data frames of primary references for the
# supported equation sets.  Exported functions: hd_citations(),
# biomass_citations(), volume_citations().

#' Citations for the supported height-diameter equation forms
#'
#' @description
#' Returns a data frame giving the primary published source for each
#' height-diameter model form supported by \code{\link{ht_from_dbh}} and
#' \code{\link{fit_hd_model}}.
#'
#' @param short logical. If \code{TRUE}, returns a compact two-column table
#'   (\code{model}, \code{author}) suitable for joining to a tree-level
#'   dataset via the \code{model} column.  Default \code{FALSE}.
#'
#' @return A data frame. When \code{short = FALSE} (default): full table with
#'   columns \code{model}, \code{n_params}, \code{equation}, \code{author},
#'   \code{citation}.  When \code{short = TRUE}: two columns only,
#'   \code{model} and \code{author}.
#'
#' @export
#' @examples
#' hd_citations()
#' hd_citations(short = TRUE)  # compact join table
hd_citations <- function(short = FALSE) {
  out <- data.frame(
    model = c(
      "naslund", "curtis", "logistic",
      "korf", "weibull", "richards"
    ),
    n_params = c(2L, 2L, 3L, 3L, 3L, 3L),
    equation = c(
      "H = 1.3 + DBH^2 / (a + b*DBH)^2",
      "H = 1.3 + a * (DBH / (1 + DBH))^b",
      "H = 1.3 + a / (1 + b*exp(-c*DBH))",
      "H = 1.3 + a * exp(-b * DBH^(-c))",
      "H = 1.3 + a * (1 - exp(-b * DBH^c))",
      "H = 1.3 + a * (1 - exp(-b*DBH))^c"
    ),
    author = c(
      "Naslund (1937)",
      "Curtis (1967)",
      "logistic (various)",
      "Korf (1939)",
      "Weibull / Yang et al. (1978)",
      "Richards (1959)"
    ),
    citation = c(
      paste0(
        "Naslund, M. (1937). Skogsforsokanstaltens gallringsfoersok i tallskog. ",
        "Meddelanden fran Statens Skogsfoersoeksanstalt, 29, 1-169."
      ),
      paste0(
        "Curtis, R.O. (1967). Height-diameter and height-diameter-age equations for ",
        "second-growth Douglas-fir. Forest Science, 13(4), 365-375."
      ),
      paste0(
        "Standard three-parameter logistic function; no single primary source. ",
        "Commonly used in forestry H-D modelling (e.g. Huang et al. 1992, ",
        "Canadian Journal of Forest Research, 22, 1146-1158)."
      ),
      paste0(
        "Korf, V. (1939). Prispevek k matematicke definici vzrustoveho zakona ",
        "lesnich porostu. Lesnicka Prace, 18, 339-356."
      ),
      paste0(
        "Yang, R.C., Kozak, A., Smith, J.H.G. (1978). The potential of Weibull-type ",
        "functions as flexible growth curves. Canadian Journal of Forest Research, ",
        "8(4), 424-431."
      ),
      paste0(
        "Richards, F.J. (1959). A flexible growth function for empirical use. ",
        "Journal of Experimental Botany, 10(29), 290-300."
      )
    ),
    stringsAsFactors = FALSE
  )
  if (short) out[, c("model", "author")] else out
}


#' Citations for the supported biomass equation sources
#'
#' @description
#' Returns a data frame giving the full bibliographic reference for each
#' \code{paper_source} value supported by \code{\link{biomass_tree}} and
#' \code{\link{biomass_components}}.
#'
#' @param short logical. If \code{TRUE}, returns a compact two-column table
#'   (\code{paper_source}, \code{author}) suitable for joining to a tree-level
#'   dataset via the \code{paper_source} column.  Default \code{FALSE}.
#'
#' @return A data frame. When \code{short = FALSE} (default): full table with
#'   columns \code{paper_source}, \code{author}, \code{species_covered},
#'   \code{citation}.  When \code{short = TRUE}: two columns only,
#'   \code{paper_source} and \code{author}.
#'
#' @export
#' @examples
#' biomass_citations()
#' biomass_citations(short = TRUE)  # compact join table
biomass_citations <- function(short = FALSE) {
  # Count species per source from the package data, minus "generic"
  counts <- tapply(
    biomass_coefs$species[biomass_coefs$species != "generic"],
    biomass_coefs$paper_source[biomass_coefs$species != "generic"],
    function(x) length(unique(x))
  )

  out <- data.frame(
    paper_source = c("Lambert2005", "Ung2008"),
    author = c(
      "Lambert, Ung & Raulier (2005)",
      "Ung, Bernier & Guo (2008)"
    ),
    species_covered = as.integer(c(
      counts[["Lambert2005"]],
      counts[["Ung2008"]]
    )),
    citation = c(
      paste0(
        "Lambert, M.-C., Ung, C.-H., Raulier, F. (2005). Canadian national tree ",
        "aboveground biomass equations. Canadian Journal of Forest Research, ",
        "35(8), 1996-2018. https://doi.org/10.1139/x05-112"
      ),
      paste0(
        "Ung, C.-H., Bernier, P., Guo, X.-J. (2008). Canadian national biomass ",
        "equations: new parameter estimates that include British Columbia data. ",
        "Canadian Journal of Forest Research, 38(5), 1123-1132. ",
        "https://doi.org/10.1139/X07-224"
      )
    ),
    stringsAsFactors = FALSE
  )
  if (short) out[, c("paper_source", "author")] else out
}


#' Citation for the volume taper equation source
#'
#' @description
#' Returns a data frame giving the primary published sources for the Kozak
#' taper equations implemented in \code{\link{tree_volume}}.
#'
#' @param short logical. If \code{TRUE}, returns a compact two-column table
#'   (\code{taper_eq}, \code{author}) suitable for joining to a tree-level
#'   dataset via the \code{taper_eq} column (values \code{"KBEC"} or
#'   \code{"KFIZ3"}, matching the \code{taper_eq} argument of
#'   \code{\link{tree_volume}}).  Default \code{FALSE}.
#'
#' @return A data frame. When \code{short = FALSE} (default): full table with
#'   columns \code{taper_eq}, \code{component}, \code{author},
#'   \code{citation}.  When \code{short = TRUE}: two columns only,
#'   \code{taper_eq} and \code{author}.
#'
#' @export
#' @examples
#' volume_citations()
#' volume_citations(short = TRUE)  # compact join table
#'
#' # Join to a tree dataset
#' trees <- data.frame(DBH = 25, WSV = 0.4, taper_eq = "KBEC")
#' merge(trees, volume_citations(short = TRUE), by = "taper_eq")
volume_citations <- function(short = FALSE) {
  out <- data.frame(
    taper_eq  = c("KBEC", "KFIZ3", NA_character_),
    component = c(
      "Taper equation form (KBEC) -- Kozak 2002/2004",
      "Taper equation form (KFIZ3) -- Kozak sqrt-transform form",
      "BC coefficient tables and implementation (both KBEC and KFIZ3)"
    ),
    author = c(
      "Kozak (2004)",
      "Kozak (1988)",
      "Luo / BC MoF (FAIBBase)"
    ),
    citation = c(
      paste0(
        "Kozak, A. (2004). My last words on taper equations. ",
        "The Forestry Chronicle, 80(4), 507-515. ",
        "https://doi.org/10.5558/tfc80507-4. ",
        "The published form of the KBEC equation (cube-root transform; ",
        "P = 1.3/H fixed inflection). Coefficients used in BCallometryR ",
        "originate from the 2002 BC MoF internal report ",
        "(_CBEC_2002 / _ERB_2002 arrays in vol_setup SAS macro)."
      ),
      paste0(
        "Kozak, A. (1988). A variable-exponent taper equation. ",
        "Canadian Journal of Forest Research, 18(11), 1363-1368. ",
        "https://doi.org/10.1139/x88-213. ",
        "The earlier sqrt-transform form (X = (1-sqrt(Z))/(1-sqrt(P)), ",
        "P fitted per species/FIZ zone) used in the KFIZ3 equations. ",
        "KFIZ3 coefficients are from the BC MoF vol_setup SAS macro ",
        "(A_FIZ array). Use KFIZ3 for pre-BEC BC inventory data (FRI ",
        "data with FIZ zone codes, roughly pre-2000)."
      ),
      paste0(
        "Luo, Y. FAIBBase R package (bcgov/FAIBBase), R/DIB_ICalculator.R. ",
        "Implements both the KBEC 2002 and KFIZ3 taper equations from the ",
        "BC Ministry of Forests volume compilation system. ",
        "https://github.com/bcgov/FAIBBase"
      )
    ),
    stringsAsFactors = FALSE
  )
  if (short) out[!is.na(out$taper_eq), c("taper_eq", "author")] else out
}
