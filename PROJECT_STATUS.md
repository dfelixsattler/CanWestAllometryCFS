# PROJECT STATUS — CanWestAllometryCFS

_Last reviewed: 2026-09-29_

An R package for growth-and-yield modellers, inventory foresters, and forestry
researchers to **impute tree heights**, **estimate stem volume**, and
**estimate above-ground biomass** from field data (PSP / non-PSP style),
provided the input matches the required column format.

---

## 1. Objectives (reference)

1. **Usable by end users** (G&Y modellers, inventory foresters, researchers) to
   impute heights, estimate volume, and estimate above-ground biomass from
   field data that fits a required format.
2. **Methods sourced from operational compilers** — BC MoF PSP/non-PSP
   compilation routines (FAIBCompiler / FAIBBase) **and** Alberta routines
   (volume functions adopted from the GYPSY model).
3. **Designed to integrate into workflows** — modular, scriptable functions
   with a consistent column contract.

---

## 2. What is done (general)

**Package health**
- Renamed `BCallometryRCFS` → `CanWestAllometryCFS`.
- `R CMD check`: 0 errors / 0 warnings / 0 notes.
- 233 unit tests passing (`testthat` edition 3).
- 4 vignettes build cleanly; README with pipeline diagram + column reference.
- Citation infrastructure (`inst/CITATION`, `citation()`), Apache 2.0 + Crown
  copyright.

**BC pipeline (Objective 2 — BC side: complete)**
- Species crosswalk: `species_correction()` → `bc_species_to_sp0()` /
  `bc_species_to_sp_type()` / `bc_species_to_biomass_name()`.
- H-D modelling: `fit_hd_models_by_group()` (SP0 → SP_TYPE fallback,
  mixed/fixed decision) + `ht_impute()` (BLUP-calibrated, broken-top option,
  multi-year QC). Validated against FAIBCompiler.
- Volume: `tree_volume()` (Kozak KBEC 16 species × 13 BEC zones; KFIZ3 for
  pre-BEC data), `tree_profile()`, `tree_volume_section()`. WSV verified to
  machine precision vs FAIBBase.
- Broken-top recovery: `ht_from_btop()` (taper-based).

**Biomass (national — works Canada-wide)**
- `biomass_tree()` / `biomass_components()`: Lambert (2005) + Ung (2008)
  national equations, 45 species, per-component and per-row height/DBH switch.

**Alberta / GYPSY (Objective 2 — AB side: started)**
- `tree_volume_huang()`: Huang (1994) variable-exponent taper volume
  (merch + whole-stem), stratified by species × Alberta natural subregion,
  configurable utilization standards (default 13/7/30). Verified to machine
  precision against the source `V_myHuang`.
- `taper_coefs_huang` dataset (8 species × subregions).

---

## 3. Outstanding items (most → least important)

### P1 — Blocking the "western Canada" claim / core objectives

1. **Integrate Alberta volume into the workflow (Objective 3).**
   `tree_volume_huang()` currently requires `genus.species` codes
   (e.g. `"PINU.CON"`) supplied directly and an Alberta natural-subregion code.
   There is **no crosswalk** from raw inventory codes to Huang codes and **no
   subregion helper**. Needed: an `ab_species_*` (or unified) crosswalk and a
   lat/long-or-code → natural-subregion lookup so AB data flows through the
   same pattern as BC.

2. **Alberta height imputation (Objective 1) is missing.**
   Only *volume* was adopted from GYPSY; there is no AB H-D model / height
   imputation. Users must supply heights for AB trees. Decide whether to (a)
   port GYPSY's H-D model, or (b) document that `fit_hd_models_by_group()` can
   be trained on AB data directly (it is region-agnostic) and provide an AB
   example.

3. **Update `DESCRIPTION` Description + package docs for western Canada.**
   The `Description:` field and `@description` still say "British Columbia
   forestry workflows" and "16 BC species"; the title was changed but the body
   was not. Reader-facing inconsistency.

### P2 — Discoverability & usability

4. **Alberta is invisible in README and vignettes.**
   No mention of Alberta / Huang / GYPSY in README; no Alberta vignette. Add
   (a) a README section + pipeline note, and (b) a `volume_ab_workflow.Rmd`
   vignette mirroring the BC volume vignette.

5. **Document the required input column contract in one place.**
   The column reference exists in README; promote it to a short "data format"
   vignette or a helper (`validate_input()`) so users can confirm their data
   fits before running — directly supports Objective 1 ("as long as input fits
   the required format").

6. **SK / MB guidance.** The CIPHA workflow applied Huang (Alberta subregions)
   to SK/MB plots. Document the recommended subregion choice for SK/MB, or add
   an explicit mapping, so those provinces are supported intentionally rather
   than by improvisation.

### P3 — Breadth & polish

7. **Expand Huang species coverage.** Only 8 species are present
   (`ABIE.BAL, PICE.GLA, PICE.MAR, PINU.CON, POPU.BAL, POPU.TRE, UNKN.HWD,
   UNKN.SWD`). Confirm whether the full Huang (1994) / GYPSY species set should
   be included; wire in the softwood/hardwood `UNKN.*` fallback logic.

8. **Unified species code strategy.** BC uses `SP0` (e.g. `"F"`, `"S"`); AB
   uses `genus.species`. Consider a single crosswalk that accepts raw codes
   from either jurisdiction and routes to the correct volume backend.

9. **GitHub repo + local folder rename.** GitHub repo is still
   `BCallometryRCFS` (needs Settings → Rename to `CanWestAllometryCFS`; then
   `git remote set-url`). Local folder is still `C:\BCallometryR`. README /
   CITATION URLs already point to the new name.

10. **Pre-existing biomass scope note.** Biomass is national (works AB/SK/MB),
    but there is no AB-specific species-name helper equivalent to
    `bc_species_to_biomass_name()`. Add an AB → common-name mapping for a
    seamless AB biomass path.

11. **Add tests for the Alberta workflow end-to-end** once the crosswalk /
    subregion helpers exist (currently only unit-level Huang volume tests).

---

## 4. Where the package deviates from its stated goals

| Goal | Status | Deviation |
|------|--------|-----------|
| Impute heights from field data | BC done / AB not yet | No AB H-D model; AB heights must be supplied. |
| Estimate volume from field data | BC done / AB partial | AB volume works but is not reachable via a crosswalk; needs species + subregion helpers. |
| Estimate above-ground biomass | done (national) | No AB-specific species-name helper; must map manually. |
| Methods from BC + Alberta compilers | BC done / AB partial | Only GYPSY *volume* adopted; other GYPSY components (H-D) not ported. |
| Integrate into workflows | BC done / AB not yet | AB feature is standalone; does not yet fit the crosswalk -> model -> estimate pattern. |
| Reflect western-Canada scope in docs | not done | Title updated; `Description`, `@description`, README, vignettes still BC-framed. |

---

## 5. Suggested next sprint (smallest set that moves the needle)

1. Update `DESCRIPTION` `Description:` + package `@description` for western
   Canada (P1-3) — quick, removes the biggest inconsistency.
2. Add an Alberta species/subregion crosswalk + `volume_ab_workflow.Rmd`
   vignette (P1-1, P2-4) — makes the AB feature actually usable in a workflow.
3. Document AB height handling: show `fit_hd_models_by_group()` trained on AB
   data (P1-2) — closes the height-imputation gap without new modelling.
4. Rename the GitHub repo + local folder (P3-9) — housekeeping.
