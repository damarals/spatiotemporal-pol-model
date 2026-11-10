# TROPOMI design: quadrature, supports, the standardized contrast, rows and the
# collapse into pseudo-rows.

column_case <- function(unit = "new-traffic-1", K = 4L, width_m = 8000) {
  x <- test_inputs()
  part <- partition(x, unit)
  calendar <- .months("2020-01", max(c(part$train$month, part$held$month, part$column$month)))
  basis <- profile_basis(x$acag_cells, K, width_m)
  list(x = x, part = part, calendar = calendar, basis = basis,
    design = column_design(part$column, x$column_cells, x$acag_cells, calendar, basis))
}

testthat::test_that("the quadrature integrates each cell exactly", {
  skip_without_moqa()
  x <- test_inputs()
  q <- .column_quadrature(x$column_cells)
  testthat::expect_identical(length(q$cell_ids), 22L)
  testthat::expect_identical(nrow(q$points), 22L * 9L)
  testthat::expect_equal(unname(drop(rowsum(q$weight, q$cell))), rep(1, 22L), tolerance = 1e-12)
  # A bilinear cell's area is its shoelace area; the order-3 rule is exact for it.
  shoelace <- vapply(q$cell_ids, function(id) {
    c4 <- x$column_cells[x$column_cells$cell_id == id, ]
    c4 <- c4[match(c("SW", "SE", "NE", "NW"), c4$corner), ]
    xs <- c4$x_m - c4$x_m[1]
    ys <- c4$y_m - c4$y_m[1]
    0.5 * abs(sum(xs * c(ys[-1], ys[1]) - c(xs[-1], xs[1]) * ys))
  }, numeric(1))
  testthat::expect_equal(q$area_m2, unname(shoelace), tolerance = 1e-10)
  # Each cell's points lie inside its corners' bounding box.
  for (i in seq_along(q$cell_ids)) {
    c4 <- x$column_cells[x$column_cells$cell_id == q$cell_ids[i], ]
    pts <- q$points[q$cell == i, , drop = FALSE]
    testthat::expect_true(all(pts[, 1] > min(c4$x_m) & pts[, 1] < max(c4$x_m)))
    testthat::expect_true(all(pts[, 2] > min(c4$y_m) & pts[, 2] < max(c4$y_m)))
  }
  rule <- gauss_legendre_3
  testthat::expect_equal(rule$nodes, c(-sqrt(3 / 5), 0, sqrt(3 / 5)), tolerance = 1e-15)
  testthat::expect_equal(rule$weights, c(5, 8, 5) / 9, tolerance = 1e-15)
})

testthat::test_that("the design has the declared supports, contrasts, months and rows", {
  skip_without_moqa()
  set.seed(1)
  cc <- column_case()
  des <- cc$design
  acag <- cc$x$acag_cells
  testthat::expect_identical(c(des$n_cells, des$K, des$n_col, des$T), c(22L, 4L, 1584L, 72L))
  testthat::expect_identical(dim(des$B_col), c(22L, 3L))
  # H0(h_S) = 0 and H0(h_S^2) = 1 over the municipal cells.
  W <- .column_support(.xy(acag), cc$basis$knots, 8000)
  wbar <- drop(crossprod(acag$weight, W))
  testthat::expect_equal(wbar, cc$basis$wbar, tolerance = 1e-12)
  B <- sweep(W, 2L, wbar) %*% cc$basis$U
  B <- B / sqrt(sum(acag$weight * rowSums(B^2)))
  testthat::expect_equal(drop(crossprod(acag$weight, B)), rep(0, 3L), tolerance = 1e-12)
  testthat::expect_equal(sum(acag$weight * rowSums(B^2)), 1, tolerance = 1e-12)
  # The cell mean of h_S is linear in z_N through B_col.
  z <- stats::rnorm(3L)
  testthat::expect_equal(drop(des$B_col %*% z), dense_h_cell(list(input = cc$x,
    prep = list(basis = cc$basis), model = list(width_m = 8000)), z), tolerance = 1e-10)
  testthat::expect_gt(diff(range(des$B_col %*% z)), 0.1)
  # Rows in month order on the calendar; cell contrasts sum to zero.
  testthat::expect_false(is.unsorted(des$rows$t))
  testthat::expect_identical(des$rows$t, sort(match(cc$part$column$month, cc$calendar)))
  testthat::expect_identical(dim(des$U_cell), c(22L, 21L))
  testthat::expect_equal(colSums(des$U_cell), rep(0, 21L), tolerance = 1e-12)
  testthat::expect_equal(crossprod(des$U_cell), diag(21L), tolerance = 1e-12)
  # Eight knots with 4 km kernels.
  d8 <- column_case(K = 8L, width_m = 4000)$design
  testthat::expect_identical(dim(d8$B_col), c(22L, 7L))
  # Whole years keep only the months before the year.
  for (year in 2023:2025) {
    cy <- column_case(paste0("year-", year))
    testthat::expect_true(all(cc$calendar[cy$design$rows$t] < paste0(year, "-01")))
    testthat::expect_identical(cy$design$n_col, 22L * 12L * (year - 2020L))
  }
})

testthat::test_that("conditioning on the pseudo-rows equals conditioning on every TROPOMI row", {
  skip_without_moqa()
  set.seed(20261110)
  cc <- column_case("omitted-months")
  des <- cc$design
  pseudo <- column_pseudo_rows(des)
  testthat::expect_identical(pseudo$t, 1:72)
  rows <- des$rows
  w <- 1 / rows$d^2
  testthat::expect_equal(pseudo$ybar, unname(drop(rowsum(w * rows$y, rows$t)) / drop(rowsum(w, rows$t))),
    tolerance = 1e-12)
  # A latent AR(1) signal e_t loads every row of month t by l; the rows share q_t
  # and cell effects. Conditioning a held signal month on all rows or on the
  # pseudo-rows gives the same law.
  l <- 0.7; s <- 0.4; rho <- 0.8; sq <- 0.3; rq <- 0.5; kT <- exp(1.5)
  b <- stats::rnorm(21L, 0, 0.2)
  mu_rows <- 1.2 + drop(des$U_cell %*% b)[rows$cell]
  lag <- abs(outer(rows$t, rows$t, "-"))
  V_rows <- (l^2 * s^2 * rho^lag + sq^2 * rq^lag) + diag(rows$d^2 / kT)
  held_t <- 30L
  c_rows <- l * s^2 * rho^abs(rows$t - held_t)
  full_mean <- drop(crossprod(c_rows, solve(V_rows, rows$y - mu_rows)))
  full_var <- s^2 - drop(crossprod(c_rows, solve(V_rows, c_rows)))
  mu_p <- 1.2 + drop(pseudo$Gbar %*% b)
  lag_p <- abs(outer(pseudo$t, pseudo$t, "-"))
  V_p <- (l^2 * s^2 * rho^lag_p + sq^2 * rq^lag_p) + diag(1 / (kT * pseudo$S))
  c_p <- l * s^2 * rho^abs(pseudo$t - held_t)
  testthat::expect_equal(drop(crossprod(c_p, solve(V_p, pseudo$ybar - mu_p))), full_mean, tolerance = 1e-9)
  testthat::expect_equal(s^2 - drop(crossprod(c_p, solve(V_p, c_p))), full_var, tolerance = 1e-9)
})
