# Roxygen2 documentation for the package's lazy-loaded datasets:
# biomass_coefs, bc_species_codes, plant_codes,
# taper_coefs_kbec, taper_coefs_kfiz3, psp_trees.

#' Simulated PSP tree list for vignette examples
#'
#' @description
#' A simulated permanent-sample-plot (PSP) tree list that mirrors the
#' structure of real BC inventory data. Designed to support all three
#' CanWestAllometryCFS workflows without any additional data preparation:
#' \itemize{
#'   \item \strong{H-D modelling} -- fit \code{\link{fit_hd_model}} to the
#'     trees with measured heights, then predict heights for unmeasured trees.
#'   \item \strong{Biomass estimation} -- pass raw species codes and BEC zones
#'     to \code{\link{bc_species_to_biomass_name}}, then call
#'     \code{\link{biomass_tree}}.
#'   \item \strong{Volume estimation} -- pass species codes to
#'     \code{\link{bc_species_to_sp0}}, then call \code{\link{tree_volume}}
#'     for trees with known heights.
#' }
#'
#' Heights are measured for approximately 25\% of trees; the remainder have
#' \code{HEIGHT = NA}, consistent with real PSP measurement protocols where
#' only a subsample of trees is height-measured per visit.
#' Approximately 5\% of trees are flagged as broken-top (\code{BTOP = TRUE});
#' these must be excluded from H-D model fitting because their total height is
#' unknown. Approximately 12\% of trees are standing dead (\code{LV_D = "D"});
#' these are excluded from H-D model training but retained for volume and
#' biomass estimation.
#' \code{HEIGHT_TRUE} records the simulated true total height for every tree and
#' is useful for validating height-prediction models.
#'
#' @format A data frame with 588 rows and 9 variables:
#' \describe{
#'   \item{SITE_IDENTIFIER}{character. Plot identifier (e.g. \code{"PSP-001"}).
#'     20 plots (5 per BEC zone).}
#'   \item{BEC_ZONE}{character. BEC zone of the plot: one of \code{"SBS"},
#'     \code{"ICH"}, \code{"CWH"}, \code{"IDF"}.}
#'   \item{unitreeid}{character. Unique tree identifier
#'     (e.g. \code{"PSP-001-01"}).}
#'   \item{SPECIES}{character. Raw BC inventory species code as recorded in
#'     PSP/VRI data: one of \code{"AT"}, \code{"FDI"}, \code{"HW"},
#'     \code{"PLI"}, \code{"SW"}.  Use \code{\link{bc_species_to_sp0}} or
#'     \code{\link{bc_species_to_biomass_name}} to translate to the formats
#'     required by the allometric functions.}
#'   \item{DBH}{numeric. Diameter at breast height (cm).}
#'   \item{HEIGHT}{numeric. Measured total height (m). \code{NA} for
#'     approximately 75\% of trees (unmeasured) and for all broken-top trees.}
#'   \item{HEIGHT_TRUE}{numeric. True simulated total height (m) for all trees.
#'     Available for validation; would not exist in real inventory data.}
#'   \item{BTOP}{logical. \code{TRUE} if the tree has a broken top. Broken-top
#'     trees have no measurable total height and are excluded from H-D model
#'     fitting. Approximately 5\% of trees.}
#'   \item{LV_D}{character. Tree status: \code{"L"} = live, \code{"D"} = dead
#'     (standing snag). Dead trees are excluded from H-D model training but
#'     retained for volume and biomass estimation. Approximately 12\% of trees.}
#' }
#' @source Simulated with \code{set.seed(42)} using Naslund H-D curves with
#'   plot-level random effects; see \code{data-raw/psp_trees.R}.
#' @seealso \code{\link{bc_species_to_sp0}},
#'   \code{\link{bc_species_to_biomass_name}}, \code{\link{fit_hd_model}},
#'   \code{\link{biomass_tree}}, \code{\link{tree_volume}}
#' @examples
#' # Overview
#' str(psp_trees)
#' table(psp_trees$SPECIES)
#' table(psp_trees$BEC_ZONE)
#'
#' # Trees with measured heights
#' measured <- psp_trees[!is.na(psp_trees$HEIGHT), ]
#' nrow(measured)  # ~150
"psp_trees"

#' Aboveground biomass coefficients (Lambert 2005 / Ung 2008)
#'
#' @description
#' A data frame of allometric power-equation coefficients for computing
#' component-level aboveground tree biomass. Coefficients are transcribed from
#' Lambert et al. (2005) and Ung et al. (2008), covering 45 Canadian tree
#' species (including all major BC commercial species).
#'
#' The data are used internally by \code{\link{biomass_tree}} and
#' \code{\link{biomass_components}}, but can also be queried directly.
#'
#' @format A data frame with 412 rows and 7 variables:
#' \describe{
#'   \item{species}{character. Lowercase common species name
#'     (e.g. \code{"lodgepole pine"}, \code{"western redcedar"}).}
#'   \item{paper_source}{character. Equation source: \code{"Lambert2005"} or
#'     \code{"Ung2008"}.}
#'   \item{height_included}{logical. \code{FALSE} for DBH-only equations;
#'     \code{TRUE} for DBH + height equations.}
#'   \item{component}{character. Biomass component: \code{"wood"},
#'     \code{"bark"}, \code{"branches"}, or \code{"foliage"}.}
#'   \item{a1}{numeric. Scaling coefficient.}
#'   \item{a2}{numeric. DBH exponent.}
#'   \item{a3}{numeric. Height exponent (\code{NA} when
#'     \code{height_included = FALSE}; \code{0} when height has no effect on
#'     that component).}
#' }
#'
#' @details
#' \strong{Equations:}
#' \itemize{
#'   \item DBH-only:     \eqn{B = a_1 \cdot DBH^{a_2}}
#'   \item DBH + height: \eqn{B = a_1 \cdot DBH^{a_2} \cdot H^{a_3}}
#' }
#' Total aboveground biomass is the sum across all four components.
#'
#' \strong{Coverage:} Most species have \code{"Lambert2005"} coefficients only.
#' The following have both Lambert2005 and Ung2008 variants: black spruce,
#' lodgepole pine, trembling aspen, white birch, white spruce, hardwood
#' (generic), softwood (generic).
#'
#' The special entry \code{species = "generic"} provides fallback coefficients
#' (height-included only) for unrecognised species names.
#'
#' @references
#'   Lambert M-C, Ung C-H, Raulier F (2005). Canadian national tree aboveground
#'   biomass equations. \emph{Canadian Journal of Forest Research}, 35(8),
#'   1996–2018.
#'
#'   Ung C-H, Bernier P, Guo X-J (2008). Canadian national biomass equations:
#'   new parameter estimates that include British Columbia data.
#'   \emph{Canadian Journal of Forest Research}, 38(5), 1123–1132.
#'
#' @seealso \code{\link{biomass_tree}}, \code{\link{biomass_components}}
#' @examples
#' # Which species are covered?
#' sort(unique(biomass_coefs$species))
#'
#' # Inspect lodgepole pine coefficients
#' biomass_coefs[biomass_coefs$species == "lodgepole pine", ]
"biomass_coefs"


#' Kozak (2002) KBEC taper equation coefficients
#'
#' @description
#' A data frame of Kozak (2002) variable-exponent taper equation coefficients
#' for 16 BC tree species across 13 BEC zones (208 rows). Coefficients are
#' transcribed from the BC Ministry of Forests volume compilation system
#' (KBEC 2002 form, \code{_CBEC_2002} array in the original SAS macro), as
#' reproduced in the FAIBBase R package (\code{DIB_ICalculator.R}).
#'
#' The data are used internally by \code{\link{tree_volume}}, but can also
#' be queried directly.
#'
#' @format A data frame with 208 rows and 12 variables:
#' \describe{
#'   \item{species}{character. BC species code (SP0): one of \code{"AC"},
#'     \code{"AT"}, \code{"B"}, \code{"C"}, \code{"D"}, \code{"E"},
#'     \code{"F"}, \code{"H"}, \code{"L"}, \code{"MB"}, \code{"PA"},
#'     \code{"PL"}, \code{"PW"}, \code{"PY"}, \code{"S"}, \code{"Y"}.}
#'   \item{bec_zone}{character. BEC zone code: one of \code{"AT"},
#'     \code{"BWBS"}, \code{"CDF"}, \code{"CWH"}, \code{"ESSF"},
#'     \code{"ICH"}, \code{"IDF"}, \code{"MH"}, \code{"MS"}, \code{"PP"},
#'     \code{"SBPS"}, \code{"SBS"}, \code{"SWB"}.}
#'   \item{B1, B2, B3}{numeric. Scale, DBH-exponent, and height-exponent
#'     parameters.}
#'   \item{B4, B5, B6, B7, B8, B9}{numeric. Exponent-term coefficients.}
#'   \item{ERR}{numeric. Multiplicative bias-correction factor
#'     (\code{_ERB_2002}).}
#' }
#'
#' @details
#' The Kozak (2002) KBEC variable-exponent taper equation is:
#' \deqn{d_i = ERR \cdot B_1 \cdot D^{B_2} \cdot H^{B_3} \cdot
#'   X^{(B_4 Z^4 + B_5 e^{-D/H} + B_6 X^{0.1} + B_7/D + B_8 H^Q + B_9 X)}}
#' where \eqn{Z = h_i/H}, \eqn{P = 1.3/H}, \eqn{Q = 1 - Z^{1/3}},
#' \eqn{X = Q / (1 - P^{1/3})}.
#'
#' Multiple (species, BEC-zone) combinations share the same equation number
#' and thus identical coefficients; the dataset has already resolved those
#' look-ups so each row is ready to use directly.
#'
#' @references
#'   Kozak, A. (2004). My last words on taper equations. \emph{The Forestry
#'   Chronicle}, 80(4), 507--515. \doi{10.5558/tfc80507-4}
#'
#' @seealso \code{\link{tree_volume}}
#' @examples
#' # Which BEC zones are covered for western hemlock?
#' taper_coefs_kbec[taper_coefs_kbec$species == "H", c("bec_zone","B1","B2")]
"taper_coefs_kbec"

#' Kozak KFIZ3 taper equation coefficients (FIZ-stratified)
#'
#' @description
#' A data frame of Kozak variable-exponent (sqrt-transform form) taper
#' equation coefficients for 16 BC tree species across 12 Forest Inventory
#' Zone (FIZ) codes (192 rows).  Coefficients are transcribed from the BC
#' Ministry of Forests volume compilation system (KFIZ3 form, \code{A_FIZ}
#' array in the original SAS macro), as reproduced in the FAIBBase R package
#' (\code{DIB_ICalculator.R}).
#'
#' The data are used internally by \code{\link{tree_volume}} when
#' \code{taper_eq = "KFIZ3"}, but can also be queried directly.
#'
#' @format A data frame with 192 rows and 11 variables:
#' \describe{
#'   \item{species}{character. BC species code (SP0): one of \code{"AC"},
#'     \code{"AT"}, \code{"B"}, \code{"C"}, \code{"D"}, \code{"E"},
#'     \code{"F"}, \code{"H"}, \code{"L"}, \code{"MB"}, \code{"PA"},
#'     \code{"PL"}, \code{"PW"}, \code{"PY"}, \code{"S"}, \code{"Y"}.}
#'   \item{fiz_zone}{character. FIZ zone code: one of \code{"A"} through
#'     \code{"L"} (12 zones; \code{"I"} is not used as a FIZ zone in BC).
#'     Grouped into three equation sets: A/B/C (northern/interior),
#'     D--J (central/southern), K/L (coastal).}
#'   \item{P}{numeric. Fitted inflection parameter (0.20, 0.25, or 0.30);
#'     defines the point of maximum taper rate. Unlike KBEC where the
#'     inflection point is fixed at \eqn{P = 1.3/H}, here \code{P} is a
#'     species/zone-specific fitted value.}
#'   \item{A0}{numeric. Scale coefficient.}
#'   \item{A1}{numeric. DBH exponent.}
#'   \item{A2}{numeric. DBH-in-exponent coefficient
#'     (\eqn{A_2^{\text{DBH}}} term; values near 1 yield small adjustments).}
#'   \item{A3, A4, A5, A6, A7}{numeric. Exponent-term coefficients
#'     (see Details).}
#' }
#'
#' @details
#' The KFIZ3 variable-exponent taper equation (Kozak sqrt-transform form):
#' \deqn{d_i = A_0 \cdot D^{A_1} \cdot A_2^D \cdot X^{E_X}}
#' where
#' \deqn{X = \frac{1 - \sqrt{Z}}{1 - \sqrt{P}}, \quad Z = h_i/H,}
#' \deqn{E_X = A_3 Z^2 + A_4 \ln(Z+0.001) + A_5 \sqrt{X}
#'            + A_6 (D/H) + A_7 e^Z.}
#'
#' The shape variable \eqn{X} uses a square-root transform of relative
#' height, with the inflection point \eqn{P} fitted per species/zone.
#' This contrasts with the KBEC equation, which uses a cube-root transform
#' and a fixed inflection at \eqn{P = 1.3/H}.
#'
#' Multiple (species, FIZ-zone) combinations share the same equation number
#' within each of the three FIZ zone groups; the coefficients have already
#' been expanded so each row is ready to use directly.
#'
#' \strong{When to use KFIZ3 vs KBEC:}
#' \itemize{
#'   \item Use \strong{KBEC} (\code{\link{taper_coefs_kbec}}) for any
#'     modern BC inventory data (VRI, current PSP protocols) where BEC zone
#'     is known.  KBEC is the current standard in BC volume compilations.
#'   \item Use \strong{KFIZ3} when working with older BC Forest Resources
#'     Inventory (FRI) data (roughly pre-2000) that carries FIZ zone codes
#'     but not BEC zone codes, or when replicating older MoF volume
#'     compilations that used the FIZ system.
#' }
#'
#' @references
#'   Luo, Y. (FAIBBase R package, bcgov/FAIBBase),
#'   \code{R/DIB_ICalculator.R}.
#'   Implements the BC MoF KFIZ3 taper equations from the original SAS
#'   \code{vol_setup} macro.
#'   \url{https://github.com/bcgov/FAIBBase}
#'
#' @seealso \code{\link{taper_coefs_kbec}}, \code{\link{tree_volume}}
#' @examples
#' # Which coefficients does Douglas-fir get in FIZ zone D?
#' taper_coefs_kfiz3[taper_coefs_kfiz3$species == "F" &
#'                   taper_coefs_kfiz3$fiz_zone == "D", ]
#'
#' # All FIZ zones covered for western hemlock
#' taper_coefs_kfiz3[taper_coefs_kfiz3$species == "H",
#'                   c("fiz_zone", "P", "A0", "A1")]
"taper_coefs_kfiz3"


#' Huang (1994) Alberta taper equation coefficients
#'
#' @description
#' Coefficients for the Huang et al. (1994) ecologically based variable-exponent
#' taper equation, used by \code{\link{tree_volume_huang}} to compute
#' individual-tree merchantable and whole-stem volume for major western Canadian
#' tree species.  Coefficients are stratified by species and Alberta natural
#' subregion.
#'
#' @format A data frame with 129 rows and 11 variables:
#' \describe{
#'   \item{species}{character. Huang species code in genus.species form
#'     (e.g. \code{"PINU.CON"}, \code{"POPU.TRE"}, \code{"PICE.GLA"}).}
#'   \item{NaturalSubregionCode}{character. Alberta natural subregion code
#'     (e.g. \code{"LF"}, \code{"CM"}, \code{"Province"}).}
#'   \item{NaturalSubregionNum}{integer. Numeric natural subregion identifier.}
#'   \item{a0, a1, a2}{numeric. Multiplicative taper coefficients.}
#'   \item{b1, b2, b3, b4, b5}{numeric. Variable-exponent taper coefficients.}
#' }
#'
#' @details
#' The Huang taper equation (Kozak sqrt-transform form) predicts inside-bark
#' diameter \eqn{d_i} at height \eqn{h_i}:
#' \deqn{d_i = a_0 \cdot D^{a_1} \cdot a_2^{D} \cdot
#'   X^{(b_1 Z^2 + b_2 \ln(Z + 0.001) + b_3 \sqrt{Z} + b_4 e^{Z} + b_5 D/H)}}
#' where \eqn{Z = h_i / H} and \eqn{X = (1 - \sqrt{Z})/(1 - \sqrt{0.225})}.
#'
#' Coefficients were transcribed from the CTAE R package
#' (\code{ptompalski/CTAE}, object \code{parameters_HuangV}), which digitises
#' Huang et al. (1994) Appendix 3.  The Upper Boreal Highlands (\code{"UB"})
#' share the Lower Boreal Highlands (\code{"LBH"}) coefficients.
#'
#' @references
#'   Huang, S. (1994). Ecologically based individual tree volume estimation
#'   for major Alberta tree species. Alberta Environmental Protection, Land
#'   and Forest Services, Forest Management Division.
#'
#' @seealso \code{\link{tree_volume_huang}}, \code{\link{taper_coefs_kbec}}
#' @examples
#' # Species covered
#' unique(taper_coefs_huang$species)
#'
#' # Coefficients for lodgepole pine in the Lower Foothills subregion
#' taper_coefs_huang[taper_coefs_huang$species == "PINU.CON" &
#'                   taper_coefs_huang$NaturalSubregionCode == "LF", ]
"taper_coefs_huang"


#' Species crosswalk: OSM PlantCodes ↔ BC SP0 ↔ common names
#'
#' @description
#' A cross-reference table mapping OSM \code{PlantCodes} strings (USDA-style
#' codes used in the \href{https://github.com/OSM-Contributors/OSM}{OSM.Allometry}
#' .NET library) to BC SP0 codes and the lowercase common names used by
#' \code{\link{biomass_tree}}.
#'
#' The table covers all 16 SP0 species in \code{\link{taper_coefs_kbec}} and
#' all 45 species in \code{\link{biomass_coefs}}, including entries that
#' appear in only one of the two equation sets.  Use the crosswalk helper
#' functions (\code{\link{bc_species_to_osm}}, \code{\link{bc_species_to_sp0}},
#' \code{\link{bc_species_to_biomass_name}}) for programmatic translation.
#'
#' @section Composite BC SP0 codes:
#' Two BC SP0 codes cover multiple species:
#' \describe{
#'   \item{\code{"S"} (spruce)}{Maps to white spruce (\code{"PIGL"}) as the
#'     primary interior BC species; Sitka spruce (\code{"PISI"}) and Engelmann
#'     spruce (\code{"PIEN"}) are also listed with \code{sp0 = "S"}.}
#'   \item{\code{"AC"} (poplar/cottonwood)}{Maps to black cottonwood
#'     (\code{"POBAT"}) as the primary coastal BC species; balsam poplar
#'     (\code{"POBA2"}) is also listed with \code{sp0 = "AC"}.}
#' }
#'
#' @format A data frame with 49 rows and 4 variables:
#' \describe{
#'   \item{plant_code}{character. OSM \code{PlantCodes} enum string
#'     (e.g. \code{"TSHE"}, \code{"PICO"}).}
#'   \item{latin_name}{character. Scientific name.}
#'   \item{common_name}{character. Lowercase common name matching
#'     \code{\link{biomass_coefs}\$species}; \code{NA} for species without
#'     Lambert/Ung biomass equations.}
#'   \item{sp0}{character. BC SP0 code matching
#'     \code{\link{taper_coefs_kbec}\$species}; \code{NA} for species without
#'     a BC taper equation.}
#' }
#'
#' @source
#' OSM PlantCodes enum: C. Hennigar, NRCan (OSM.CommonModels, USDA PLANTS
#' codes).  BC SP0 codes: BC Ministry of Forests Timber Supply Review
#' conventions.  Common names: Lambert et al. (2005) and Ung et al. (2008).
#'
#' @seealso \code{\link{bc_species_to_osm}}, \code{\link{bc_species_to_sp0}},
#'   \code{\link{bc_species_codes}},
#'   \code{\link{taper_coefs_kbec}}, \code{\link{biomass_coefs}}
#' @examples
#' # Show all species that have both a taper equation and a biomass equation
#' plant_codes[!is.na(plant_codes$sp0) & !is.na(plant_codes$common_name), ]
"plant_codes"


#' Detailed BC inventory species codes with SP0 group mapping
#'
#' @description
#' A lookup table of 132 detailed BC forest inventory species codes and their
#' corresponding SP0 group codes used in the BC volume compilation system
#' (and in \code{\link{tree_volume}}).
#'
#' PSP and VRI records typically carry detailed species codes (e.g.
#' \code{"FDI"}, \code{"FDC"}, \code{"SW"}, \code{"SE"}) that must be
#' collapsed to SP0 codes before calling \code{\link{tree_volume}}.  Use
#' \code{\link{bc_species_to_sp0}} for that conversion.
#'
#' @format A data frame with 132 rows and 5 variables:
#' \describe{
#'   \item{species_code}{character. Detailed BC species code (e.g.
#'     \code{"FDI"}, \code{"SW"}, \code{"HW"}, \code{"PLI"}).}
#'   \item{sp0}{character. SP0 group code used in volume compilation
#'     (e.g. \code{"F"}, \code{"S"}, \code{"H"}, \code{"PL"}).}
#'   \item{description}{character. Common species name or group description.}
#'   \item{sp_sindex}{character. Species code used in BC site index
#'     calculations (connects to the SIndexR / SIndexBC package).}
#'   \item{sp_type}{character. \code{"C"} (conifer) or \code{"D"} (deciduous
#'     broadleaf).}
#' }
#'
#' @section SP0 group key:
#' The 16 SP0 codes and their primary species:
#' \tabular{ll}{
#'   \strong{SP0} \tab \strong{Primary species} \cr
#'   AC \tab Poplar / cottonwood (Populus spp.) \cr
#'   AT \tab Trembling aspen (Populus tremuloides) \cr
#'   B  \tab True firs (Abies spp.: amabilis, subalpine, balsam) \cr
#'   C  \tab Western redcedar (Thuja plicata) \cr
#'   D  \tab Red alder (Alnus rubra) \cr
#'   E  \tab Paper birch (Betula papyrifera) \cr
#'   F  \tab Douglas-fir (Pseudotsuga menziesii) \cr
#'   H  \tab Western hemlock (Tsuga heterophylla) \cr
#'   L  \tab Western larch (Larix occidentalis) \cr
#'   MB \tab Bigleaf maple (Acer macrophyllum) \cr
#'   PA \tab Whitebark pine (Pinus albicaulis) \cr
#'   PL \tab Lodgepole pine (Pinus contorta) \cr
#'   PW \tab Western white pine (Pinus monticola) \cr
#'   PY \tab Ponderosa pine (Pinus ponderosa) \cr
#'   S  \tab Spruce (Picea spp.: white, Sitka, Engelmann) \cr
#'   Y  \tab Yellow-cedar (Chamaecyparis nootkatensis) \cr
#' }
#'
#' @source
#' BC Ministry of Forests species coding conventions; reproduced from the
#' \code{lookup_species()} function in the FAIBBase R package
#' (bcgov/FAIBBase, Y. Luo, NRCan/BC MoF).
#'
#' @seealso \code{\link{bc_species_to_sp0}}, \code{\link{plant_codes}},
#'   \code{\link{tree_volume}}
#' @examples
#' # Collapse a PSP species column to SP0 before computing volume
#' bc_species_to_sp0(c("FDI", "FDC", "HW", "BL", "PLI", "SW", "SE"))
#'
#' # All species codes in the spruce group
#' bc_species_codes[bc_species_codes$sp0 == "S", c("species_code", "description")]
"bc_species_codes"
