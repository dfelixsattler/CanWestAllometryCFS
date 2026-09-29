test_that("tree_volume WSV is positive and reasonable", {
  v <- tree_volume(bec_zone = "CWH", species = "H", dbh = 30.7, height = 27.4)
  expect_length(v, 1)
  expect_true(v > 0.5 && v < 5)   # roughly 1-2 m³ for a 30 cm, 27 m tree
})

test_that("tree_volume scalar bec_zone/species recycle over dbh vector", {
  vols <- tree_volume("SBS", "PL",
                      dbh    = c(10, 20, 30, 40),
                      height = c( 8, 15, 22, 28))
  expect_length(vols, 4)
  expect_true(all(vols > 0))
  # Volume increases monotonically with DBH/height
  expect_true(all(diff(vols) > 0))
})

test_that("tree_volume WSV >= MER >= STUMP for same tree", {
  args <- list(bec_zone = "IDF", species = "D", dbh = 35, height = 30)
  v_wsv   <- do.call(tree_volume, c(args, list(volume_type = "WSV")))
  v_mer   <- do.call(tree_volume, c(args, list(volume_type = "MER")))
  v_stump <- do.call(tree_volume, c(args, list(volume_type = "STUMP")))
  expect_true(v_wsv >= v_mer)
  expect_true(v_wsv >= v_stump)
  expect_true(v_stump > 0)
})

test_that("tree_volume MER < WSV (stump + top always excluded)", {
  v_wsv <- tree_volume("CWH", "S", 40, 35, volume_type = "WSV")
  v_mer <- tree_volume("CWH", "S", 40, 35, volume_type = "MER")
  expect_true(v_wsv > v_mer)
})

test_that("tree_volume broken top reduces WSV", {
  v_whole  <- tree_volume("CWH", "H", 30.7, 27.4)
  v_broken <- tree_volume("CWH", "H", 30.7, 27.4, btop_height = 15)
  expect_true(v_broken < v_whole)
})

test_that("tree_volume multiple trees with vector inputs", {
  vols <- tree_volume(
    bec_zone = c("CWH", "CWH", "IDF"),
    species  = c("H",   "S",   "D"),
    dbh      = c(30.7,  42.3,  25.0),
    height   = c(27.4,  37.3,  22.0)
  )
  expect_length(vols, 3)
  expect_true(all(vols > 0))
})

test_that("tree_volume returns NA with warning for unknown species/BEC combo", {
  expect_warning(
    v <- tree_volume("ZZZ", "H", 30, 25),
    "coefficients"
  )
  expect_true(is.na(v))
})

test_that("tree_volume returns NA with warning for invalid DBH/height", {
  expect_warning(v <- tree_volume("CWH", "H", 0.5, 25), "invalid")
  expect_true(is.na(v))
  expect_warning(v <- tree_volume("CWH", "H", 30, 1.0), "invalid")
  expect_true(is.na(v))
})

test_that("taper_coefs_kbec dataset is correct structure", {
  expect_s3_class(taper_coefs_kbec, "data.frame")
  expect_equal(nrow(taper_coefs_kbec), 208L)
  expect_true(all(c("species", "bec_zone",
                    "B1","B2","B3","B4","B5","B6","B7","B8","B9","ERR") %in%
                    names(taper_coefs_kbec)))
  expect_equal(length(unique(taper_coefs_kbec$species)),  16L)
  expect_equal(length(unique(taper_coefs_kbec$bec_zone)), 13L)
  expect_true(all(taper_coefs_kbec$ERR > 1 & taper_coefs_kbec$ERR < 1.02))
})

test_that("volume_citations returns correct structure", {
  cit <- volume_citations()
  expect_s3_class(cit, "data.frame")
  expect_equal(nrow(cit), 3L)
  expect_named(cit, c("taper_eq", "component", "author", "citation"))

  # short form
  cit_s <- volume_citations(short = TRUE)
  expect_named(cit_s, c("taper_eq", "author"))
  expect_equal(nrow(cit_s), 2L)  # NA row dropped
})

# --- tree_profile tests -------------------------------------------------------

test_that("tree_profile returns a data.frame with expected columns", {
  prof <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
  expect_s3_class(prof, "data.frame")
  expect_named(prof, c("height", "dib", "vol_slice", "cumvol", "comment"))
})

test_that("tree_profile height grid starts at 0 and ends at tree height", {
  prof <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
  expect_equal(prof$height[1], 0)
  expect_equal(tail(prof$height, 1), 27.4)
})

test_that("tree_profile includes breast height (1.3 m) in grid", {
  prof <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
  expect_true(any(abs(prof$height - 1.3) < 1e-6))
})

test_that("tree_profile comment column has expected annotations", {
  prof <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
  expect_true(any(grepl("ground",              prof$comment)))
  expect_true(any(grepl("stump height",        prof$comment)))
  expect_true(any(grepl("breast height",       prof$comment)))
  expect_true(any(grepl("max merchantable",    prof$comment)))
  expect_true(any(grepl("tip",                 prof$comment)))
})

test_that("tree_profile cumvol is non-decreasing and matches tree_volume WSV", {
  prof <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
  expect_true(all(diff(prof$cumvol) >= -1e-10))   # non-decreasing
  wsv <- tree_volume("CWH", "H", 30.7, 27.4)
  # profile total volume should closely match tree_volume WSV
  expect_equal(tail(prof$cumvol, 1), wsv, tolerance = 1e-6)
})

test_that("tree_profile dib is non-negative and largest near base", {
  prof <- tree_profile("IDF", "D", dbh = 40, height = 32)
  expect_true(all(prof$dib >= 0))
  # Butt swell means max DIB should occur in the lower stem, not at breast height
  expect_true(which.max(prof$dib) <= which(abs(prof$height - 1.3) < 1e-6))
})

test_that("tree_profile broken top ends at btop_height", {
  prof_full  <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4)
  prof_btop  <- tree_profile("CWH", "H", dbh = 30.7, height = 27.4, btop_height = 15)
  expect_equal(tail(prof_btop$height, 1), 15)
  expect_true(tail(prof_btop$cumvol, 1) < tail(prof_full$cumvol, 1))
  expect_true(any(grepl("break height", prof_btop$comment)))
  expect_false(any(grepl("\\btip\\b",   prof_btop$comment)))
})

test_that("tree_profile volume segments sum correctly", {
  prof <- tree_profile("SBS", "PL", dbh = 25, height = 20)
  # Sum of vol_slice should equal total cumvol (last row)
  expect_equal(sum(prof$vol_slice), tail(prof$cumvol, 1), tolerance = 1e-10)
})

test_that("tree_profile returns NULL with warning for unknown BEC/species", {
  expect_warning(r <- tree_profile("ZZZ", "H", 30, 25), "coefficients")
  expect_null(r)
})

test_that("tree_profile returns NULL with warning for invalid inputs", {
  expect_warning(r <- tree_profile("CWH", "H", 0.5, 25), "invalid DBH")
  expect_null(r)
  expect_warning(r <- tree_profile("CWH", "H", 30, 1.0), "invalid height")
  expect_null(r)
})

test_that("tree_profile errors on vector dbh/height", {
  expect_error(tree_profile("CWH", "H", c(10, 20), c(10, 15)), "single tree")
})

# --- FAIBBase regression values (extracted 2026-06-30) ----------------------

test_that("tree_volume WSV matches FAIBBase reference values", {
  # Reference values from FAIBBase treeProfile() 2026-06-30.
  # Intact tree: no broken top, so BCallometryR and FAIBBase use identical
  # slice boundaries -> should match to < 1e-4 m3.
  expect_equal(
    tree_volume("CWH", "H", 30.7, 27.4),
    0.85515978,
    tolerance = 1e-4
  )
  # Broken-top tree: CanWestAllometryCFS intentionally >= FAIBBase by <= 0.02 m3.
  # CanWestAllometryCFS includes the last Smalian slice up to btop_height;
  # FAIBBase excludes it (HT_I_next < BTOPHeight strict inequality).
  ref_btop <- 0.71841045  # FAIBBase VOL_WSV for btop = 15 m
  bc_btop  <- tree_volume("CWH", "H", 30.7, 27.4, btop_height = 15)
  expect_true(bc_btop >= ref_btop)
  expect_true(bc_btop - ref_btop < 0.02)
})