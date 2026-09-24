suppressPackageStartupMessages(library(data.table))
devtools::load_all("C:/BCallometryR", quiet = TRUE)
set.seed(42)

trees <- fread("my_tests/faib_tree_detail_subtest.csv")
plots <- fread("my_tests/faib_plot_header.csv")
trees <- merge(trees, plots[, .(SITE_IDENTIFIER, BEC_ZONE)],
               by = "SITE_IDENTIFIER", all.x = TRUE)
trees[, DBH    := as.numeric(DBH)]
trees[, HEIGHT := as.numeric(HEIGHT)]
trees[, BTOP   := BROKEN_TOP_IND == "Y"]
trees[, SPECIES_CORR    := species_correction(SPECIES, BEC_ZONE)]
trees[, SPECIES_SP0     := bc_species_to_sp0(SPECIES_CORR)]
trees[, SPECIES_SP_TYPE := bc_species_to_sp_type(SPECIES_CORR)]

m <- trees[HEIGHT_SOURCE == "Field measured" & !BTOP &
           LV_D == "L" &
           !is.na(HEIGHT) & !is.na(DBH) & !is.na(SPECIES_SP0) &
           DBH > 0 & HEIGHT > 1.3]
m[, row_id := .I]

# Plot-level 80/20 split: all trees at a plot go to train OR test together
all_sites <- unique(m[SPECIES_SP0 == "H"][["SITE_IDENTIFIER"]])
set.seed(42)
test_sites <- sample(all_sites, size = ceiling(length(all_sites) * 0.20))
train <- m[!SITE_IDENTIFIER %in% test_sites]
test  <- m[ SITE_IDENTIFIER %in% test_sites]

hw_tr <- train[SPECIES == "HW"]
hw_te <- test[ SPECIES == "HW"]

cat("=== HW Training set ===\n")
cat("Trees:", nrow(hw_tr), "\n")
cat("Plots:", length(unique(hw_tr[["SITE_IDENTIFIER"]])), "\n")

cat("\n=== HW Holdout set ===\n")
cat("Trees:", nrow(hw_te), "\n")
cat("Plots:", length(unique(hw_te[["SITE_IDENTIFIER"]])), "\n")
cat("Holdout plots also in training:", 
    sum(unique(hw_te[["SITE_IDENTIFIER"]]) %in% unique(hw_tr[["SITE_IDENTIFIER"]])), "\n")
cat("Holdout plots NOT  in training:",
    sum(!unique(hw_te[["SITE_IDENTIFIER"]]) %in% unique(hw_tr[["SITE_IDENTIFIER"]])), "\n")

cat("\n=== Fitting H-D models ===\n")
hd <- fit_hd_models_by_group(as.data.frame(train))
cat("H (SP0) model:\n")
print(hd$summary[hd$summary$Group == "H", ], row.names = FALSE)

cat("\n=== HD_SOURCE for HW holdout predictions ===\n")
tb <- as.data.frame(hw_te)
tb[["HEIGHT"]] <- NA_real_
pr <- as.data.table(ht_impute(tb, hd))
print(table(pr[["HD_SOURCE"]], useNA = "ifany"))

# -----------------------------------------------------------------------
# Bias / RMSE on HOLDOUT
# -----------------------------------------------------------------------
pr[, HEIGHT_obs := hw_te[["HEIGHT"]]]
pr[, resid := HT_PROJ - HEIGHT_obs]
cat("\n=== HW Holdout: obs vs predicted ===\n")
cat(sprintf("bias (predicted - observed): %+.2f m  |  MAE: %.2f m  |  RMSE: %.2f m\n",
            mean(pr$resid, na.rm=TRUE),
            mean(abs(pr$resid), na.rm=TRUE),
            sqrt(mean(pr$resid^2, na.rm=TRUE))))
cat(sprintf("DBH range: %.1f -- %.1f cm  |  Height range: %.1f -- %.1f m\n",
            min(hw_te$DBH), max(hw_te$DBH),
            min(hw_te$HEIGHT), max(hw_te$HEIGHT)))
cat(sprintf("Predictions with BLUP: %d  |  Fixed-effects only (BLUP=0): %d\n",
            sum(!is.na(pr$HT_PROJ) & pr$SITE_IDENTIFIER %in% train$SITE_IDENTIFIER),
            sum(!is.na(pr$HT_PROJ) & !pr$SITE_IDENTIFIER %in% train$SITE_IDENTIFIER)))

# -----------------------------------------------------------------------
# Bias / RMSE on TRAINING (in-sample fit)
# -----------------------------------------------------------------------
cat("\n=== HW Training: in-sample fit (observed vs fitted) ===\n")
tr_blank <- as.data.frame(hw_tr)
tr_blank[["HEIGHT"]] <- NA_real_
pr_tr <- as.data.table(ht_impute(tr_blank, hd))
pr_tr[, HEIGHT_obs := hw_tr[["HEIGHT"]]]
pr_tr[, resid := HT_PROJ - HEIGHT_obs]
cat(sprintf("bias: %+.2f m  |  MAE: %.2f m  |  RMSE: %.2f m\n",
            mean(pr_tr$resid, na.rm=TRUE),
            mean(abs(pr_tr$resid), na.rm=TRUE),
            sqrt(mean(pr_tr$resid^2, na.rm=TRUE))))

