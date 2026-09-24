# Species crosswalk helpers ---------------------------------------------------
# Translates between detailed inventory species codes, group-level species
# codes (SP0 in BC inventory), OSM PlantCodes, and the lowercase common names
# used in biomass_coefs.
#
# Terminology used throughout:
#   detailed inventory code -- species-level code as recorded in PSP/VRI data
#                              (e.g. "FDI", "FDC", "PLI", "SW", "SE", "HW")
#   group-level code (SP0)  -- broader grouping used by the volume equations;
#                              collapses varieties into a single group
#                              (e.g. FDI + FDC -> "F"; SW + SE + SS -> "S")
#   PlantCode               -- four-to-six-character botanical code used in the
#                              OSM.Allometry library (e.g. "TSHE", "PICO")
#   common name             -- lowercase English common name as used in
#                              biomass_coefs (e.g. "western hemlock")
#
# The underlying lookup table is the plant_codes dataset.  All functions
# return NA_character_ (with a warning) for unrecognised inputs.

# Internal lookup -------------------------------------------------------------

.pc_lookup <- function(from_col, to_col, value, prefer_primary = TRUE) {
  # Returns the first match in plant_codes[from_col] == value → to_col.
  # For composite SP0 codes (e.g. S → PIGL/PISI/PIEN) the first row in the
  # table is the primary/preferred species; prefer_primary keeps that order.
  idx <- which(plant_codes[[from_col]] == value)
  if (length(idx) == 0L) return(NA_character_)
  if (prefer_primary) idx <- idx[1L]
  plant_codes[[to_col]][idx]
}


# Exported helpers ------------------------------------------------------------

#' Convert group-level species codes (SP0) to OSM PlantCodes
#'
#' @description
#' Translates one or more group-level species codes to OSM
#' \code{PlantCodes} strings (as used in the \code{OSM.Allometry} .NET
#' library).  Group-level codes collapse multiple detailed inventory codes
#' into a single group; in BC inventory this system is called SP0
#' (e.g. \code{"F"} for all Douglas-fir variants, \code{"S"} for all spruce).
#' These are the codes required by \code{\link{tree_volume}} and
#' \code{\link{taper_coefs_kbec}}.
#'
#' Composite group codes that cover multiple species (e.g. \code{"S"} for all
#' spruce, \code{"AC"} for cottonwood/poplar) return the primary species listed
#' in \code{\link{plant_codes}}: white spruce (\code{"PIGL"}) for \code{"S"},
#' black cottonwood (\code{"POBAT"}) for \code{"AC"}.
#'
#' @param sp0 character vector of group-level species codes (SP0).
#'
#' @return character vector of OSM \code{PlantCodes} strings, same length as
#'   \code{sp0}. \code{NA} is returned with a warning for unrecognised codes.
#'
#' @seealso \code{\link{plant_codes}}, \code{\link{plant_code_to_sp0}},
#'   \code{\link{common_name_to_plant_code}}
#' @noRd
.sp0_to_plant_code <- function(sp0) {
  result <- vapply(sp0, .pc_lookup,
                   from_col = "sp0", to_col = "plant_code",
                   FUN.VALUE = character(1L))
  missing <- is.na(result)
  if (any(missing)) {
    warning(".sp0_to_plant_code: unrecognised group-level code(s): ",
            paste(sp0[missing], collapse = ", "))
  }
  unname(result)
}


#' Convert OSM PlantCodes to group-level species codes (SP0)
#'
#' @description
#' Translates one or more OSM \code{PlantCodes} strings to group-level
#' species codes (referred to as SP0 in BC inventory) used in
#' \code{\link{tree_volume}}.  Returns \code{NA} for species
#' that do not have a BC taper equation in \code{\link{taper_coefs_kbec}}.
#'
#' @param plant_code character vector of OSM \code{PlantCodes} strings.
#'
#' @return character vector of group-level species codes (SP0), same length as
#'   \code{plant_code}.  \code{NA} is returned (with a warning) for
#'   unrecognised codes or species with no BC taper equation.
#'
#' @seealso \code{\link{plant_codes}}, \code{\link{bc_species_to_osm}}
plant_code_to_sp0 <- function(plant_code) {
  result <- vapply(plant_code, .pc_lookup,
                   from_col = "plant_code", to_col = "sp0",
                   FUN.VALUE = character(1L))
  missing <- is.na(result)
  if (any(missing)) {
    warning("plant_code_to_sp0: no group-level (SP0) code for: ",
            paste(plant_code[missing], collapse = ", "),
            " (species may not have a BC taper equation)")
  }
  unname(result)
}


#' Convert biomass common names to OSM PlantCodes
#'
#' @description
#' Translates one or more lowercase common species names (as used in
#' \code{\link{biomass_tree}} and \code{\link{biomass_coefs}}) to OSM
#' \code{PlantCodes} strings.
#'
#' @param common_name character vector of lowercase common names.
#'
#' @return character vector of OSM \code{PlantCodes} strings, same length as
#'   \code{common_name}. \code{NA} is returned with a warning for
#'   unrecognised names.
#'
#' @seealso \code{\link{plant_codes}}, \code{\link{bc_species_to_osm}}
common_name_to_plant_code <- function(common_name) {
  result <- vapply(tolower(trimws(common_name)), .pc_lookup,
                   from_col = "common_name", to_col = "plant_code",
                   FUN.VALUE = character(1L))
  missing <- is.na(result)
  if (any(missing)) {
    warning("common_name_to_plant_code: unrecognised common name(s): ",
            paste(common_name[missing], collapse = ", "))
  }
  unname(result)
}


#' Convert OSM PlantCodes to biomass common names
#'
#' @description
#' Translates one or more OSM \code{PlantCodes} strings to the lowercase
#' common names used in \code{\link{biomass_tree}}.  Returns \code{NA} for
#' species that do not have biomass equations in this package.
#'
#' @param plant_code character vector of OSM \code{PlantCodes} strings.
#'
#' @return character vector of common names, same length as
#'   \code{plant_code}. \code{NA} is returned (with a warning) for
#'   unrecognised codes or species with no biomass equation.
#'
#' @seealso \code{\link{plant_codes}}, \code{\link{bc_species_to_osm}}
plant_code_to_common_name <- function(plant_code) {
  result <- vapply(plant_code, .pc_lookup,
                   from_col = "plant_code", to_col = "common_name",
                   FUN.VALUE = character(1L))
  missing <- is.na(result)
  if (any(missing)) {
    warning("plant_code_to_common_name: no biomass common name for: ",
            paste(plant_code[missing], collapse = ", "),
            " (species may not have a Lambert/Ung biomass equation)")
  }
  unname(result)
}


#' Convert detailed inventory species codes to group-level codes (SP0)
#'
#' @description
#' Translates one or more detailed inventory species codes (e.g. \code{"FDI"},
#' \code{"FDC"}, \code{"SW"}, \code{"SE"}, \code{"EP"}) to group-level
#' species codes used by \code{\link{tree_volume}} and
#' \code{\link{taper_coefs_kbec}}.
#'
#' Group-level codes collapse multiple detailed variants into a single group:
#' for example, \code{"FDI"} and \code{"FDC"} both map to \code{"F"}
#' (Douglas-fir); \code{"SW"}, \code{"SE"}, and \code{"SS"} all map to
#' \code{"S"} (spruce).  In BC inventory this grouping system is called SP0.
#'
#' This is the most common translation needed in PSP and VRI workflows, where
#' tree records carry detailed species codes but the volume equations work at
#' the group level.
#'
#' @param species_code character vector of detailed inventory species codes
#'   (case-sensitive; use uppercase as in BC inventory records).
#'
#' @return character vector of group-level species codes (SP0), same length as
#'   \code{species_code}.
#'   \code{NA} is returned with a warning for unrecognised codes.
#'
#' @export
#' @seealso \code{\link{bc_species_codes}}, \code{\link{bc_species_to_osm}},
#'   \code{\link{tree_volume}}
#' @examples
#' bc_species_to_sp0(c("FDI", "FDC", "SW", "SE", "SS", "HW", "BL", "PLI"))
#' # "F"   "F"   "S"   "S"   "S"   "H"   "B"   "PL"
bc_species_to_sp0 <- function(species_code) {
  result <- bc_species_codes$sp0[
    match(species_code, bc_species_codes$species_code)]
  missing <- is.na(result)
  if (any(missing)) {
    warning("bc_species_to_sp0: unrecognised species code(s): ",
            paste(species_code[missing], collapse = ", "))
  }
  result
}


#' Convert BC inventory species codes to species type (conifer / deciduous)
#'
#' @description
#' Translates one or more detailed inventory species codes to their
#' broadest functional type: \code{"C"} (conifer) or \code{"D"} (deciduous
#' broadleaf). This is the coarsest level of the FAIBCompiler H-D model
#' hierarchy and is used as a last-resort fallback when there are too few
#' measured trees to fit a species- or SP0-level model.
#'
#' @param species_code character vector of detailed inventory species codes
#'   (case-sensitive; use uppercase as in BC inventory records).
#'
#' @return character vector of species-type codes (\code{"C"} or \code{"D"}),
#'   same length as \code{species_code}.
#'   \code{NA} is returned with a warning for unrecognised codes.
#'
#' @export
#' @seealso \code{\link{bc_species_to_sp0}}, \code{\link{bc_species_codes}}
#' @examples
#' bc_species_to_sp_type(c("FDI", "PLI", "HW", "AT", "SW"))
#' # "C"   "C"   "C"   "D"   "C"
bc_species_to_sp_type <- function(species_code) {
  result <- bc_species_codes$sp_type[
    match(species_code, bc_species_codes$species_code)]
  missing <- is.na(result)
  if (any(missing)) {
    warning("bc_species_to_sp_type: unrecognised species code(s): ",
            paste(species_code[missing], collapse = ", "))
  }
  result
}


#' Convert BC inventory species codes to OSM PlantCodes
#'
#' @description
#' Translates one or more BC inventory species codes to
#' \href{https://github.com/bcgov/FAIBBase}{OSM} \code{PlantCodes} -- the
#' four-to-six-character botanical identifiers used by the
#' \code{OSM.Allometry} .NET library.
#'
#' Internally maps through the group-level code (SP0), so multiple detailed
#' codes that share a group return the same primary \code{PlantCode}
#' (e.g. \code{"FDI"} and \code{"FDC"} both return \code{"PSME"};
#' \code{"SW"}, \code{"SE"}, and \code{"SS"} all return \code{"PIGL"}
#' -- the primary spruce).
#'
#' @param species_code character vector of detailed BC inventory species codes
#'   (uppercase, as recorded in PSP or VRI data).
#'
#' @return character vector of OSM \code{PlantCodes} strings, same length as
#'   \code{species_code}.  \code{NA} is returned with a warning for
#'   unrecognised codes or species with no OSM entry.
#'
#' @export
#' @seealso \code{\link{bc_species_codes}}, \code{\link{plant_codes}},
#'   \code{\link{bc_species_to_sp0}}
#' @examples
#' bc_species_to_osm(c("FDI", "FDC", "HW", "BL", "PLI", "SW"))
bc_species_to_osm <- function(species_code) {
  sp0 <- suppressWarnings(bc_species_to_sp0(species_code))
  result <- suppressWarnings(.sp0_to_plant_code(sp0))
  missing <- is.na(result)
  if (any(missing)) {
    warning("bc_species_to_osm: no OSM PlantCode for: ",
            paste(species_code[missing], collapse = ", "))
  }
  result
}


#' Correct BC inventory species codes based on BEC zone context
#'
#' @description
#' Standardises raw detailed inventory species codes using the BEC-zone-aware
#' correction rules from the BC Ministry of Forests volume compilation system.
#'
#' Many detailed codes are ambiguous without BEC context:
#' \itemize{
#'   \item Generic \code{"S"} becomes \code{"SW"}, \code{"SE"}, or \code{"SS"}
#'         depending on zone.
#'   \item Generic \code{"B"} / \code{"BA"} becomes \code{"BL"} (interior) or
#'         \code{"BA"} / \code{"BG"} (coastal).
#'   \item Generic \code{"H"} becomes \code{"HW"} or \code{"HM"} depending on
#'         zone.
#' }
#' In addition, several shorthand codes are unconditionally remapped (e.g.
#' \code{"FDI"} / \code{"FDC"} \eqn{\to} \code{"FD"}, \code{"L"}
#' \eqn{\to} \code{"LW"}, \code{"P"} / \code{"PLI"} \eqn{\to} \code{"PL"}).
#'
#' The corrected codes can then be passed to \code{\link{bc_species_to_sp0}}
#' to obtain the group-level code (SP0) needed by \code{\link{tree_volume}}.
#'
#' @param species character vector of raw detailed inventory species codes
#'   (uppercase, as recorded in PSP or VRI data).
#' @param bec_zone character vector of BEC zone codes (e.g. \code{"CWH"},
#'   \code{"SBS"}, \code{"ESSF"}).  Recycled to \code{length(species)} if
#'   scalar.
#' @param bec_subzone character vector of BEC subzone codes (e.g.
#'   \code{"mm"}, \code{"wk"}).  Only required for the \code{CWHmm}
#'   edge-case (B/BL \eqn{\to} BG); defaults to \code{""}.  Recycled to
#'   \code{length(species)} if scalar.
#'
#' @return character vector of corrected species codes, same length as
#'   \code{species}.
#'
#' @details
#' \strong{Correction order (mirrors the BC MoF species code correction routine):}
#' \enumerate{
#'   \item Unconditional one-to-one remaps (e.g. \code{FDI/FDC} \eqn{\to}
#'         \code{FD}, \code{L} \eqn{\to} \code{LW}, \code{HX} \eqn{\to}
#'         \code{HW}).
#'   \item Codes still longer than 2 characters are truncated to 2
#'         (e.g. \code{HXM} \eqn{\to} \code{HX}).
#'   \item BEC-based disambiguation (spruce, balsam, hemlock).
#' }
#' Note: \code{HX} in the \emph{input} is always mapped to \code{HW} by rule 1
#' (before BEC logic), so it will not be further corrected to \code{HM} even
#' in MH or ESSF zones.  \code{HXM} goes through truncation first, then BEC
#' disambiguation, so it correctly resolves to \code{HM} in those zones.
#'
#' @references
#'   Luo, Y. (BC MoF / NRCan). \emph{FAIBCompiler R package},
#'   \code{speciesCorrection()} function. bcgov/FAIBCompiler (GitHub).
#'   Correction rules attributed to Rene and Dan, email May 12, 2021.
#'
#' @export
#' @seealso \code{\link{bc_species_to_sp0}}, \code{\link{bc_species_codes}},
#'   \code{\link{tree_volume}}
#' @examples
#' # Typical PSP/VRI pipeline:
#' species_raw <- c("FDI", "FDC", "HXM", "B", "S", "SXE", "PLI", "L")
#' bec         <- c("CWH", "IDF", "MH",  "SBS","CWH","ESSF","SBS","ICH")
#' corrected   <- species_correction(species_raw, bec)
#' corrected   # "FD"  "FD"  "HM"  "BL" "SS"  "SE"  "PL"  "LW"
#' bc_species_to_sp0(corrected)
#' #  "F"   "F"   "H"   "B"  "S"   "S"   "PL"  "L"
species_correction <- function(species, bec_zone, bec_subzone = "") {
  n <- length(species)
  if (length(bec_zone)    == 1L) bec_zone    <- rep(bec_zone,    n)
  if (length(bec_subzone) == 1L) bec_subzone <- rep(bec_subzone, n)

  sp <- species

  # 1. Unconditional one-to-one remaps (order matches FAIBCompiler) ----------
  sp[sp == "A"]                      <- "AT"
  sp[sp == "AX"]                     <- "AC"
  sp[sp == "HX"]                     <- "HW"   # before BEC disambiguation
  sp[sp == "C"]                      <- "CW"
  sp[sp %in% c("D", "RA")]           <- "DR"
  sp[sp %in% c("E", "EXP", "EA")]    <- "EP"
  sp[sp == "J"]                      <- "JR"
  sp[sp == "L"]                      <- "LW"
  sp[sp %in% c("P", "PLI")]          <- "PL"
  sp[sp %in% c("FDI", "FDC")]        <- "FD"
  sp[sp %in% c("SX", "SXL", "SXW")] <- "SW"
  sp[sp == "SXE"]                    <- "SE"
  sp[sp == "T"]                      <- "TW"
  sp[sp %in% c("X", "XC")]           <- "XC"
  sp[sp == "ZH"]                     <- "XH"

  # 2. Truncate remaining codes > 2 characters (e.g. HXM → HX) -------------
  long <- nchar(sp) > 2L
  sp[long] <- substr(sp[long], 1L, 2L)

  # 3. BEC-based disambiguation ---------------------------------------------
  is_coastal <- bec_zone %in% c("CWH", "CDF", "MH", "CMA")
  bec_i_c    <- ifelse(is_coastal, "C", "I")
  done       <- logical(n)

  # -- Spruce --
  # S / SE in interior (not ESSF) → SW
  sel <- bec_i_c == "I" & bec_zone != "ESSF" & sp %in% c("S", "SE") & !done
  sp[sel] <- "SW"; done[sel] <- TRUE

  # S in ESSF → SE (interior Engelmann zone)
  sel <- bec_zone == "ESSF" & sp == "S" & !done
  sp[sel] <- "SE"; done[sel] <- TRUE

  # SS in interior → SW
  sel <- bec_i_c == "I" & sp == "SS" & !done
  sp[sel] <- "SW"; done[sel] <- TRUE

  # S on coast → SS (Sitka)
  sel <- bec_i_c == "C" & sp == "S" & !done
  sp[sel] <- "SS"; done[sel] <- TRUE

  # SW / SE in CWH → SS
  sel <- bec_zone == "CWH" & sp %in% c("SE", "SW") & !done
  sp[sel] <- "SS"; done[sel] <- TRUE

  # SW in ESSF / MH → SE
  sel <- bec_zone %in% c("ESSF", "MH") & sp == "SW" & !done
  sp[sel] <- "SE"; done[sel] <- TRUE

  # -- Balsam / true firs --
  # B / BA in interior → BL (subalpine fir)
  sel <- bec_i_c == "I" & sp %in% c("B", "BA") & !done
  sp[sel] <- "BL"; done[sel] <- TRUE

  # B / BL in CWHmm → BG (grand fir; not legitimate on coast otherwise)
  sel <- bec_zone == "CWH" & bec_subzone == "mm" & sp %in% c("B", "BL") & !done
  sp[sel] <- "BG"; done[sel] <- TRUE

  # B / BL on coast → BA (amabilis fir)
  sel <- bec_i_c == "C" & sp %in% c("B", "BL") & !done
  sp[sel] <- "BA"; done[sel] <- TRUE

  # BG in interior (not ICH) → BL
  sel <- bec_i_c == "I" & bec_zone != "ICH" & sp == "BG" & !done
  sp[sel] <- "BL"; done[sel] <- TRUE

  # BG on coast (not CWH / CDF) → BA
  sel <- bec_i_c == "C" & !(bec_zone %in% c("CWH", "CDF")) & sp == "BG" & !done
  sp[sel] <- "BA"; done[sel] <- TRUE

  # -- Hemlock --
  # H / HXM (→ HX after truncation) in MH / ESSF → HM (mountain hemlock)
  sel <- bec_zone %in% c("MH", "ESSF") & sp %in% c("H", "HXM", "HX") & !done
  sp[sel] <- "HM"; done[sel] <- TRUE

  # H / HXM (→ HX) on coast or in ICH / SBS → HW (western hemlock)
  sel <- (bec_i_c == "C" | bec_zone %in% c("ICH", "SBS")) &
         sp %in% c("H", "HXM", "HX") & !done
  sp[sel] <- "HW"; done[sel] <- TRUE

  # -- Misc --
  # VP on coast → VB (bitter cherry)
  sel <- bec_i_c == "C" & sp == "VP" & !done
  sp[sel] <- "VB"

  sp
}


#' Convert BC inventory species codes to biomass common names
#'
#' @description
#' Translates BC inventory species codes to the lowercase common names
#' expected by \code{\link{biomass_tree}}, applying the BEC-zone-aware
#' \code{\link{species_correction}} step automatically.
#'
#' This is the recommended entry point for biomass workflows: supply the
#' raw species codes and BEC zone exactly as they appear in inventory data
#' (PSP, VRI, or any other BC inventory source) and the function handles
#' BEC disambiguation and species-name resolution internally.
#'
#' The full conversion pipeline is:
#' \enumerate{
#'   \item \code{\link{species_correction}} — standardise and BEC-disambiguate
#'         codes (e.g. generic \code{"S"} → \code{"SW"}, \code{"SE"}, or
#'         \code{"SS"} depending on zone; \code{"B"} → \code{"BL"} or
#'         \code{"BA"}; \code{"H"} → \code{"HW"} or \code{"HM"}).
#'   \item Common-name resolution -- a small set of corrected codes would
#'         otherwise resolve to the wrong primary species via the standard
#'         group-level lookup (e.g. \code{"SE"} and \code{"SS"} share the
#'         spruce group with white spruce; \code{"BL"} shares the balsam group
#'         with amabilis fir).  These are resolved directly from the corrected
#'         detailed code.  All other codes go through
#'         \code{\link{bc_species_to_sp0}} ->
#'         internal group-level -> PlantCode lookup ->
#'         \code{\link{plant_code_to_common_name}}.
#' }
#'
#' @param species     character vector of detailed BC inventory species codes
#'   (uppercase, as recorded in the database).
#' @param bec_zone    character vector of BEC zone codes (e.g.
#'   \code{"SBS"}, \code{"CWH"}).  Recycled to \code{length(species)} if
#'   scalar.
#' @param bec_subzone character vector of BEC subzone codes.  Only needed
#'   for the \code{CWHmm} grand-fir edge-case; defaults to \code{""}.
#'
#' @return character vector of biomass common names, same length as
#'   \code{species}, suitable for passing to \code{\link{biomass_tree}}.
#'   \code{NA} is returned (with a warning) for species that have no
#'   matching entry.
#'
#' @export
#' @seealso \code{\link{species_correction}}, \code{\link{bc_species_to_sp0}},
#'   \code{\link{biomass_tree}}, \code{\link{biomass_coefs}}
#' @examples
#' bc_species_to_biomass_name(
#'   species  = c("PLI", "SW", "SE",  "BL",  "BA",   "HW",  "FDI",  "AT"),
#'   bec_zone = c("SBS","SBS","ESSF","SBS","CWH",  "CWH", "IDF",  "SBS")
#' )
#' # "lodgepole pine"  "white spruce"  "engelmann spruce"  "alpine fir"
#' # "pacific silver fir"  "western hemlock"  "douglas-fir"  "trembling aspen"
bc_species_to_biomass_name <- function(species, bec_zone, bec_subzone = "") {
  # Step 1: BEC-aware code standardisation
  corrected <- species_correction(species, bec_zone, bec_subzone)

  # Direct overrides for corrected codes whose group-level lookup maps to the
  # wrong primary species in the plant_codes table:
  #   SP0 "S" primary → white spruce (PIGL), so SE and SS need overrides.
  #   SP0 "B" primary → pacific silver fir (ABAM), so BL needs an override.
  overrides <- c(
    SE = "engelmann spruce",
    SS = "sitka spruce",
    BL = "alpine fir"        # Abies lasiocarpa; Lambert (2005) calls it "alpine fir"
  )

  result <- vapply(corrected, function(sp) {
    if (!is.na(sp) && sp %in% names(overrides)) {
      return(overrides[[sp]])
    }
    sp0 <- suppressWarnings(bc_species_to_sp0(sp))
    if (is.na(sp0)) return(NA_character_)
    pc  <- suppressWarnings(.sp0_to_plant_code(sp0))
    if (is.na(pc)) return(NA_character_)
    suppressWarnings(plant_code_to_common_name(pc))
  }, FUN.VALUE = character(1L))

  missing <- is.na(result)
  if (any(missing)) {
    warning("bc_species_to_biomass_name: no biomass name for: ",
            paste(species[missing], collapse = ", "))
  }
  unname(result)
}
