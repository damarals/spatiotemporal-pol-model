test_that("Figure 2 draws one tile per source and month with data", {
  skip_without_moqa()
  input <- test_inputs()
  tiles <- ggplot2::ggplot_build(.source_coverage_plot(input))$data[[1L]]
  expect_equal(nrow(tiles), 230L)
  # Rows from the top: TROPOMI, SatPM2.5, reference NO2, reference PM2.5, low-cost network.
  expect_equal(as.vector(table(factor(tiles$y, levels = 5:1))), c(72L, 60L, 49L, 43L, 6L))
  expect_false(anyDuplicated(tiles[c("x", "y")]) > 0)
  expect_true(all(tiles$x %in% 1:72))
  calendar <- format(seq(as.Date("2020-01-01"), as.Date("2025-12-01"), by = "month"), "%Y-%m")
  expect_equal(range(calendar[tiles$x[tiles$y == 4]]), c("2020-01", "2024-12"))
  expect_equal(range(calendar[tiles$x[tiles$y == 1]]), c("2024-08", "2025-01"))
})
