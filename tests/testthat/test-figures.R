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

test_that("Figure 4 gives the four colour-bar legends and the two pair titles one height", {
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"), width = 7.2, height = 2.9, family = paper_font)
  withr::defer(grDevices::dev.off())
  square <- function(x) sf::st_polygon(list(cbind(x + c(0, 0.01, 0.01, 0, 0), c(0, 0, 0.01, 0.01, 0))))
  cells <- expand.grid(i = 0:2, pollutant = c("NO2", "PM25"), quantity = c("mean", "relative_width"),
    stringsAsFactors = FALSE)
  squares <- sf::st_sf(cells, value = cells$i + 1, geometry = sf::st_sfc(lapply(cells$i * 0.01, square), crs = 4326))
  outline <- sf::st_sfc(sf::st_polygon(list(cbind(c(0, 0.03, 0.03, 0, 0), c(0, 0, 0.01, 0.01, 0)))), crs = 4326)
  points <- sf::st_sf(location_id = "a", geometry = sf::st_sfc(sf::st_point(c(0.015, 0.005)), crs = 4326))
  mm <- function(u) grid::convertHeight(u, "mm", valueOnly = TRUE)
  # Equal guide rows put the bars, their titles and tick labels at one height once patchwork aligns the panels.
  guide_row <- function(p, q) {
    g <- ggplot2::ggplotGrob(.map_panel(squares, outline, points, p, q))
    mm(g$heights[g$layout$t[g$layout$name == "guide-box-bottom"]])
  }
  rows <- c(guide_row("NO2", "mean"), guide_row("NO2", "relative_width"),
    guide_row("PM25", "mean"), guide_row("PM25", "relative_width"))
  expect_equal(rows, rep(rows[1L], 4L), tolerance = 1e-9)
  # Top-aligned titles share their baseline when their heights and descents agree.
  extent <- function(p) {
    text <- grid::textGrob(.pair_title(p), gp = grid::gpar(fontfamily = paper_font, fontsize = 11.5))
    c(mm(grid::grobHeight(text)), mm(grid::grobDescent(text)))
  }
  expect_equal(extent("PM25"), extent("NO2"), tolerance = 1e-9)
})
