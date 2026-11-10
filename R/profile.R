# Spatial profiles: K knots placed as maximin points among the SatPM2.5 cells,
# Gaussian kernels of width width_m normalized at each point, normalized Helmert
# contrasts of the knot levels and the municipal reference weights wbar.

#' Profile basis shared by both pollutants.
#'
#' @param acag_cells SatPM2.5 cells with x_m, y_m and area weights.
#' @param K number of knots (4, or 8 in the sensitivity analysis).
#' @param width_m kernel width in metres (8000, or 4000).
#' @return list: knots (K x 2), U (K x K-1 contrasts), wbar (municipal weights,
#'   sum one), K, width_m.
profile_basis <- function(acag_cells, K, width_m) {
  knots <- .maximin_knots(.xy(acag_cells), acag_cells$weight, K)
  wbar <- drop(crossprod(acag_cells$weight, profile_support(.xy(acag_cells), knots, width_m)))
  if (abs(sum(wbar) - 1) > 1e-10 || any(wbar <= 0))
    stop("The municipal reference weights must be positive and sum to one.", call. = FALSE)
  list(knots = knots, U = contrasts(K), wbar = wbar, K = K, width_m = width_m)
}

#' Kernel weights of points (n x 2 matrix) on the knots, normalized per point.
profile_support <- function(points, knots, width_m) {
  kernel <- exp(-(outer(points[, 1], knots[, 1], "-")^2 + outer(points[, 2], knots[, 2], "-")^2) /
    width_m^2)
  unname(kernel / rowSums(kernel))
}

#' Normalized Helmert contrasts: K x (K-1) orthonormal columns that sum to zero.
contrasts <- function(K) {
  U <- stats::contr.helmert(K)
  unname(sweep(U, 2L, sqrt(colSums(U^2)), "/"))
}

# log r_k of a normalized profile: a = U z sd, log r_k = a_k - log sum_j wbar_j exp(a_j).
.log_profile <- function(basis, z, sd) {
  a <- drop(basis$U %*% z) * sd
  a <- a - max(a)
  a - log(sum(basis$wbar * exp(a)))
}

.xy <- function(x) as.matrix(x[, c("x_m", "y_m")])

# First knot nearest the weighted centroid, then each next one farthest from those placed.
.maximin_knots <- function(points, weights, K) {
  n <- nrow(points)
  centroid <- drop(crossprod(weights, points))
  selected <- which.min(rowSums(sweep(points, 2, centroid)^2))
  nearest <- rep(Inf, n)
  for (i in seq_len(K)) {
    nearest <- pmin(nearest, rowSums(sweep(points, 2, points[selected[i], ])^2))
    nearest[selected] <- -Inf
    if (i < K) selected <- c(selected, which.max(nearest))
  }
  points[selected, , drop = FALSE]
}
