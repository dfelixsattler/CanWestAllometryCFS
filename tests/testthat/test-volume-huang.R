# Tests for the Huang (1994) Alberta taper volume function.

test_that("tree_volume_huang returns a data.frame with merch and total", {
  v <- tree_volume_huang("PINU.CON", dbh = 25, height = 20, subregion = "LF")
  expect_s3_class(v, "data.frame")
  expect_named(v, c("merch", "total"))
  expect_equal(nrow(v), 1L)
})

test_that("tree_volume_huang gives positive, ordered volumes for a sound tree", {
  v <- tree_volume_huang("PINU.CON", dbh = 25, height = 20,
                         subregion = "LF", min_dbh = 0)
  expect_true(v$merch > 0)
  expect_true(v$total > v$merch)   # total includes tip + stump
})

test_that("min_dbh nulls out volume for small trees", {
  v <- tree_volume_huang("PINU.CON", dbh = 10, height = 12,
                         subregion = "LF", min_dbh = 13)
  expect_true(is.na(v$merch))
  expect_true(is.na(v$total))
})

test_that("min_dbh = 0 allows volume for merchantable-sized small trees", {
  # DBH must exceed the top diameter for merchantable volume to exist
  v <- tree_volume_huang("PINU.CON", dbh = 15, height = 14,
                         subregion = "LF", min_dbh = 0)
  expect_true(is.finite(v$merch))
})

test_that("tree_volume_huang is vectorised and recycles scalars", {
  v <- tree_volume_huang(
    species   = c("PINU.CON", "POPU.TRE", "PICE.GLA"),
    dbh       = c(25, 30, 40),
    height    = c(20, 24, 28),
    subregion = "LF")
  expect_equal(nrow(v), 3L)
  expect_true(all(is.finite(v$total)))
})

test_that("unknown species/subregion returns NA with a warning", {
  expect_warning(
    v <- tree_volume_huang("ZZZ.ZZZ", dbh = 25, height = 20,
                           subregion = "LF", min_dbh = 0))
  expect_true(is.na(v$merch))
})

test_that("UB subregion falls back to LBH coefficients", {
  v_ub  <- tree_volume_huang("PINU.CON", 25, 20, subregion = "UB",  min_dbh = 0)
  v_lbh <- tree_volume_huang("PINU.CON", 25, 20, subregion = "LBH", min_dbh = 0)
  expect_equal(v_ub$total, v_lbh$total)
})

test_that("larger top-diameter limit reduces merchantable volume", {
  v7  <- tree_volume_huang("PINU.CON", 30, 24, subregion = "LF",
                           min_dbh = 0, utop_dib = 7)
  v13 <- tree_volume_huang("PINU.CON", 30, 24, subregion = "LF",
                           min_dbh = 0, utop_dib = 13)
  expect_true(v13$merch < v7$merch)
})
