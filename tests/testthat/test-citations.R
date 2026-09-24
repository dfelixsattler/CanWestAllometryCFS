# Tests for the citation helper functions.

test_that("hd_citations returns full and short tables", {
  full <- hd_citations()
  expect_s3_class(full, "data.frame")
  expect_true(all(c("model", "n_params", "equation", "author", "citation") %in%
                    names(full)))

  short <- hd_citations(short = TRUE)
  expect_named(short, c("model", "author"))
  expect_equal(nrow(short), nrow(full))
})

test_that("biomass_citations returns full and short tables", {
  full <- biomass_citations()
  expect_s3_class(full, "data.frame")
  expect_true(all(c("paper_source", "author", "citation") %in% names(full)))

  short <- biomass_citations(short = TRUE)
  expect_named(short, c("paper_source", "author"))
})

test_that("volume_citations returns full and short tables", {
  full <- volume_citations()
  expect_s3_class(full, "data.frame")
  expect_true("taper_eq" %in% names(full))

  short <- volume_citations(short = TRUE)
  expect_named(short, c("taper_eq", "author"))
  expect_true(all(!is.na(short$taper_eq)))
})
