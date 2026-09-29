# CanWestAllometryCFS

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/dfelixsattler/CanWestAllometryCFS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dfelixsattler/CanWestAllometryCFS/actions/workflows/R-CMD-check.yaml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
<!-- badges: end -->

An R package for tree-level calculation of H-D relationships and height
imputation, stem volume, and above-ground biomass, with an emphasis on
compatibility with BC Ministry of Forests' PSP and non-PSP field data.

Many of the package's functions were designed to replicate the BC Forest
Analysis and Inventory Branch's (FAIB) compilation routines for PSP and
non-PSP data, found in FAIBCompiler and FAIBBase. The current package
attempts to improve on the existing functions by (1) making them more
accessible through simplified functions and comprehensive vignettes,
(2) giving the user greater control over key decisions (e.g. model form
selection, broken-top tree handling, mixed- vs. fixed-effects fitting), and
(3) improving documentation.

While this package was designed with BC PSP and non-PSP data in mind, plot
data from other jurisdictions can also be used. The main requirement is that
column names in the user's dataset match those expected by the functions (or
that the relevant `_col` arguments are set to the correct column names).

## Installation

```r
# Install from GitHub
remotes::install_github("dfelixsattler/CanWestAllometryCFS")
```

## Data pipeline

The diagram below shows the end-to-end workflow from raw field data to
estimated heights, volumes, and biomass.

```mermaid
flowchart TD
    A[("Field data\nDBH · HEIGHT · SPECIES · BEC_ZONE · LV_D")]

    A --> SC

    subgraph SC["1 · Species code crosswalk"]
        direction TB
        SC1["species_correction(species, bec_zone)\ntranslate & BEC-disambiguate raw codes"]
        SC2["bc_species_to_sp0()\nSP0 group code  ·  e.g. F, H, S, PL"]
        SC3["bc_species_to_sp_type()\nSP_TYPE  ·  C (conifer) or D (deciduous)"]
        SC4["bc_species_to_biomass_name()\ncommon name  ·  e.g. douglas-fir"]
        SC1 --> SC2
        SC1 --> SC3
        SC1 --> SC4
    end

    SC2 --> HD

    subgraph HD["2 · H-D modelling  (alive, measured trees only)"]
        direction TB
        HD1["fit_hd_models_by_group(measured)\nNaslund / Curtis · mixed-effects or fixed-effects\nSP0 level → SP_TYPE fallback"]
        HD2["ht_impute(all_trees, hd_result, impute_btop = TRUE)\nHT_PROJ for every tree\nBLUP-calibrated per plot"]
        HD1 --> HD2
    end

    SC4 --> BIO
    HD2 --> VOL
    HD2 --> BIO

    subgraph VOL["3 · Volume"]
        VOL1["tree_volume(bec_zone, sp0, dbh, ht_proj)\nWSV · MER · NMR · STUMP  (m³)\nKozak KBEC or KFIZ3 taper equation"]
    end

    subgraph BIO["4 · Biomass"]
        BIO1["biomass_tree(species_name, dbh, height)\nTotal aboveground dry biomass (kg)\nLambert 2005 / Ung 2008 national equations"]
        BIO2["biomass_components()\nwood · bark · branches · foliage (kg)"]
    end
```

> **Notes**
> - Dead trees (`LV_D == "D"`) and broken-top trees (`BTOP == TRUE`) are
>   excluded from H-D model training but flow through to volume and biomass.
> - `impute_btop = TRUE` estimates total height for broken-top trees from DBH,
>   matching the BC MoF compilation routine.
> - Use `taper_eq = "KFIZ3"` with FIZ zone codes for pre-BEC inventory data.

## Quick start

```r
library(CanWestAllometryCFS)

trees <- psp_trees

# 1. Species crosswalk
trees$SPECIES_CORR    <- species_correction(trees$SPECIES, trees$BEC_ZONE)
trees$SPECIES_SP0     <- bc_species_to_sp0(trees$SPECIES_CORR)
trees$SPECIES_SP_TYPE <- bc_species_to_sp_type(trees$SPECIES_CORR)
trees$SPECIES_NAME    <- bc_species_to_biomass_name(trees$SPECIES, trees$BEC_ZONE)

# 2. Fit H-D models and impute heights
measured  <- trees[!is.na(trees$HEIGHT) & trees$LV_D == "L" & !trees$BTOP, ]
hd_result <- fit_hd_models_by_group(measured)
trees     <- ht_impute(trees, hd_result, impute_btop = TRUE)

# 3. Volume
trees$WSV_M3 <- tree_volume(trees$BEC_ZONE, trees$SPECIES_SP0,
                             trees$DBH, trees$HT_PROJ)

# 4. Biomass
trees$BIOMASS_KG <- biomass_tree(trees$SPECIES_NAME, trees$DBH,
                                  height = trees$HT_PROJ)
```

## Vignettes

| Vignette | Topic |
|---|---|
| `vignette("full_psp_workflow")` | End-to-end PSP workflow |
| `vignette("hd_psp_workflow")` | H-D modelling in depth |
| `vignette("biomass_psp_workflow")` | Biomass equations |
| `vignette("volume_psp_workflow")` | Kozak taper volume |

## Column reference

### Raw PSP data columns (inputs)

The data passed to the package functions should have, at minimum, the
following columns:

| Column | Description |
|---|---|
| `SITE_IDENTIFIER` | Plot/site identifier |
| `VISIT_NUMBER` | Measurement period (multi-year data) |
| `unitreeid` | Unique tree identifier across all visits |
| `SPECIES` | Raw BC inventory species code (e.g. `"PLI"`, `"FDI"`, `"HW"`) |
| `BEC_ZONE` | Biogeoclimatic Ecosystem Classification zone (e.g. `"SBS"`, `"CWH"`) |
| `DBH` | Diameter at breast height (cm) |
| `HEIGHT` | Measured total height (m); `NA` for unmeasured trees |
| `BTOP` | `TRUE` if the tree has a broken top |
| `LV_D` | Tree status: `"L"` = live, `"D"` = dead standing |

### Columns added by the package

| Column | Added by | Description |
|---|---|---|
| `SPECIES_CORR` | `species_correction()` | BEC-disambiguated species code |
| `SPECIES_SP0` | `bc_species_to_sp0()` | SP0 group code used by H-D and volume functions |
| `SPECIES_SP_TYPE` | `bc_species_to_sp_type()` | `"C"` (conifer) or `"D"` (deciduous); H-D fallback grouping |
| `SPECIES_NAME` | `bc_species_to_biomass_name()` | Common name used by biomass equations |
| `HT_PROJ` | `ht_impute()` | Filled height (m); measured or model-estimated |
| `HD_SOURCE` | `ht_impute()` | Source of the height value (e.g. `"measured"`, `"naslund_mixed-effects_sp0"`) |
| `HT_FLAG` | `ht_impute()` | QC flag (e.g. `"btop"` for broken-top trees) |

## How to Cite

If you use CanWestAllometryCFS in published work, please cite it as:

> Sattler, D. (2026). *CanWestAllometryCFS: Allometric Equations for British Columbia: Height–Diameter, Taper, Volume, and Biomass*. R package version 0.1.0. https://github.com/dfelixsattler/CanWestAllometryCFS

A machine-readable citation is also available in R via:

```r
citation("CanWestAllometryCFS")
```

## Acknowledgements

The allometric methods implemented in this package build on the work of
**Yong Luo** (BC Ministry of Forests) and the
[FAIBBase](https://github.com/bcgov/FAIBBase) and
[FAIBCompiler](https://github.com/bcgov/FAIBCompiler) packages.

## License

© His Majesty the King in Right of Canada, as represented by the Minister of
Natural Resources, 2026.

This package is distributed under the [Apache License 2.0](LICENSE).

Portions of this package are derived from
[FAIBBase](https://github.com/bcgov/FAIBBase) and
[FAIBCompiler](https://github.com/bcgov/FAIBCompiler),
© 2019 Province of British Columbia, which are themselves distributed under the
Apache License 2.0. Those portions remain subject to that license, and the
original copyright and attribution notices are retained here as required by
Section 4 of the Apache License 2.0.
