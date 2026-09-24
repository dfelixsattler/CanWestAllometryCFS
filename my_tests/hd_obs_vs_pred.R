# =============================================================================
# HD model validation — observed vs predicted
#   1. BCallometryR obs vs predicted (80/20 train/test on measured trees)
#   2. BCallometryR vs FAIB fixed-effects on same holdout (fair comparison)
#
# Dataset: faib_tree_detail_subtest.csv + faib_plot_header.csv
# =============================================================================
devtools::load_all("C:/BCallometryR", quiet = TRUE)
suppressPackageStartupMessages(library(data.table))

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Load and prepare
# -----------------------------------------------------------------------------
trees <- fread("my_tests/faib_tree_detail_subtest.csv")
plots <- fread("my_tests/faib_plot_header.csv")

trees <- merge(trees, plots[, .(SITE_IDENTIFIER, BEC_ZONE)],
               by = "SITE_IDENTIFIER", all.x = TRUE)
trees[, DBH    := as.numeric(DBH)]
trees[, HEIGHT := as.numeric(HEIGHT)]
trees[, BTOP   := BROKEN_TOP_IND == "Y"]
# LV_D is already in the file ("L" = live, "D" = dead)

# Species crosswalk
trees[, SPECIES_CORR    := species_correction(SPECIES, BEC_ZONE)]
trees[, SPECIES_SP0     := bc_species_to_sp0(SPECIES_CORR)]
trees[, SPECIES_SP_TYPE := bc_species_to_sp_type(SPECIES_CORR)]

# Restrict to field-measured, alive, non-btop, valid DBH/height
measured_all <- trees[HEIGHT_SOURCE == "Field measured" &
                      !BTOP & LV_D == "L" &
                      !is.na(HEIGHT) & !is.na(DBH) & !is.na(SPECIES_SP0) &
                      DBH > 0 & HEIGHT > 1.3]

cat(sprintf("Field-measured trees available: %d across %d sites\n\n",
            nrow(measured_all), length(unique(measured_all$SITE_IDENTIFIER))))

# -----------------------------------------------------------------------------
# 1. 80/20 train/test split (stratified by SP0 so rare species are covered)
# -----------------------------------------------------------------------------
measured_all[, row_id := .I]
test_idx <- measured_all[, .(row_id = sample(row_id, max(1L, ceiling(.N * 0.20)))),
                          by = SPECIES_SP0]$row_id

train <- measured_all[!row_id %in% test_idx]
test  <- measured_all[ row_id %in% test_idx]

cat(sprintf("Train: %d trees  |  Test (holdout): %d trees\n\n",
            nrow(train), nrow(test)))

# -----------------------------------------------------------------------------
# 2. Fit BCallometryRCFS H-D models on training set
# -----------------------------------------------------------------------------
cat("=== Fitting BCallometryRCFS H-D models on training set ===\n")
hd_result <- fit_hd_models_by_group(as.data.frame(train))
print(hd_result$summary, row.names = FALSE)
cat("\n")

# -----------------------------------------------------------------------------
# 3. Predict for holdout set
#    Blank the measured heights first — ht_impute() only predicts where NA
# -----------------------------------------------------------------------------
test_blank          <- as.data.frame(test)
test_blank$HEIGHT   <- NA_real_   # hide the true heights before predicting
test_pred           <- ht_impute(test_blank, hd_result)
test_pred           <- as.data.table(test_pred)

# Add actual observed heights back for comparison
test_pred[, HEIGHT_obs := test$HEIGHT]

# Only evaluate trees that got a prediction
eval <- test_pred[!is.na(HT_PROJ) & !is.na(HEIGHT_obs)]

cat(sprintf("Holdout trees with prediction: %d  |  No prediction: %d\n\n",
            nrow(eval),
            nrow(test_pred) - nrow(eval)))

# -----------------------------------------------------------------------------
# 4. BCallometryRCFS: observed (measured HEIGHT) vs predicted (HT_PROJ)
# -----------------------------------------------------------------------------
cat("=== BCallometryRCFS: observed vs predicted (holdout set) ===\n")
eval[, resid := HT_PROJ - HEIGHT_obs]

cat(sprintf("Overall  — bias: %+.2f m  |  MAE: %.2f m  |  RMSE: %.2f m\n\n",
            mean(eval$resid),
            mean(abs(eval$resid)),
            sqrt(mean(eval$resid^2))))

sp_res <- eval[, .(n       = .N,
                   HT_mean = round(mean(HEIGHT_obs), 1),
                   bias    = round(mean(resid), 2),
                   mae     = round(mean(abs(resid)), 2),
                   rmse    = round(sqrt(mean(resid^2)), 2)),
               by = .(SPECIES, SPECIES_SP0)][order(SPECIES)]
cat("By species:\n")
print(sp_res, row.names = FALSE)
cat("\n")

# -----------------------------------------------------------------------------
# 5. Compare: BCallometryRCFS vs FAIB fixed-effects on same holdout
#    FAIB fixed-effects = population-average prediction (no BLUP)
#    extracted from hd_result fixed coefficients, applied without site BLUP
# -----------------------------------------------------------------------------
cat("=== Comparison: BCallometryRCFS (with BLUP) vs fixed-effects only ===\n")
cat("Fixed-effects only = population-average prediction (no plot-level BLUP),\n")
cat("mimicking what FAIB would give for a site not in its training data.\n\n")

# Generate fixed-effects-only predictions by zeroing the BLUP
# We do this by imputing on a copy of test where SITE_IDENTIFIER is set to a
# value not in the training data (so BLUP lookup returns 0)
test_nosite                    <- as.data.frame(test)
test_nosite$HEIGHT             <- NA_real_          # blank heights
test_nosite$SITE_IDENTIFIER    <- "UNSEEN_SITE"     # no BLUP available → fixed effects only
pred_fixedonly                 <- ht_impute(test_nosite, hd_result)
pred_fixedonly                 <- as.data.table(pred_fixedonly)

eval2 <- data.table(
  SPECIES      = eval$SPECIES,
  SPECIES_SP0  = eval$SPECIES_SP0,
  HEIGHT_obs   = eval$HEIGHT_obs,
  HT_BLUP      = eval$HT_PROJ,
  HT_FIXED     = pred_fixedonly$HT_PROJ[match(eval$row_id, test$row_id)]
)
eval2 <- eval2[!is.na(HT_BLUP) & !is.na(HT_FIXED)]

eval2[, resid_blup  := HT_BLUP  - HEIGHT_obs]
eval2[, resid_fixed := HT_FIXED - HEIGHT_obs]

cat(sprintf("With BLUP     — bias: %+.2f m  |  MAE: %.2f m  |  RMSE: %.2f m\n",
            mean(eval2$resid_blup),  mean(abs(eval2$resid_blup)),
            sqrt(mean(eval2$resid_blup^2))))
cat(sprintf("Fixed-eff only— bias: %+.2f m  |  MAE: %.2f m  |  RMSE: %.2f m\n\n",
            mean(eval2$resid_fixed), mean(abs(eval2$resid_fixed)),
            sqrt(mean(eval2$resid_fixed^2))))

sp_comp <- eval2[, .(
  n          = .N,
  HT_mean    = round(mean(HEIGHT_obs), 1),
  blup_mae   = round(mean(abs(resid_blup)),  2),
  fixed_bias = round(mean(resid_fixed), 2),
  fixed_mae  = round(mean(abs(resid_fixed)), 2)
), by = .(SPECIES, SPECIES_SP0)][order(SPECIES)]
cat("By species (BLUP vs fixed-effects only):\n")
print(sp_comp, row.names = FALSE)
