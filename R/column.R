# TROPOMI column on the NO2 profile's basis (natural scale, unit 1e-5 mol/m2):
# cell averages of the standardized NO2 contrast over an order-3 Gauss-Legendre
# quadrature of each whole cell, cell contrasts, and the collapse of each
# month's cells into one precision-weighted pseudo-row.

# Order-3 Gauss-Legendre rule on [-1, 1], as the published fits computed it.
gauss_legendre_3 <- list(nodes = c(-0.77459666924148352, 4.2500725161431774e-17, 0.7745966692414834),
  weights = c(0.55555555555555558, 0.88888888888888873, 0.5555555555555558))

#' Column design of one partition.
#'
#' @param column TROPOMI rows available to the fit (partition()$column).
#' @param column_cells cell corners (read_inputs()$column_cells).
#' @param acag_cells SatPM2.5 cells, the municipal reference.
#' @param calendar model months "YYYY-MM".
#' @param basis profile_basis() of the NO2 profile.
#' @return list: rows (sorted by month and cell: cell, t, y, d), n_col, n_cells,
#'   cell_ids, K, B_col (n_cells x K-1 cell means of the standardized
#'   contrast), U_cell (cell contrasts), T.
column_design <- function(column, column_cells, acag_cells, calendar, basis) {
  t <- match(column$month, calendar)
  if (!nrow(column) || anyNA(t)) stop("Column rows must lie in the model calendar.", call. = FALSE)
  if (any(!is.finite(column$value)) || any(!(column$sd > 0)))
    stop("Column values must be finite and their SDs positive.", call. = FALSE)
  quadrature <- .column_quadrature(column_cells)
  A_col <- rowsum(.column_support(quadrature$points, basis$knots, basis$width_m) * quadrature$weight,
    quadrature$cell)
  dimnames(A_col) <- NULL
  W_acag <- .column_support(.xy(acag_cells), basis$knots, basis$width_m)
  wbar <- drop(crossprod(acag_cells$weight, W_acag))
  B_acag <- sweep(W_acag, 2L, wbar) %*% basis$U
  r_B <- sqrt(sum(acag_cells$weight * rowSums(B_acag^2)))
  cell <- match(column$cell_id, quadrature$cell_ids)
  if (anyNA(cell)) stop("A column row names an unknown cell.", call. = FALSE)
  rows <- data.frame(cell = cell, t = t, y = column$value, d = column$sd)
  rows <- rows[order(rows$t, rows$cell), , drop = FALSE]
  rownames(rows) <- NULL
  n_cells <- length(quadrature$cell_ids)
  list(rows = rows, n_col = nrow(rows), n_cells = n_cells, cell_ids = quadrature$cell_ids,
    K = basis$K, B_col = sweep(A_col, 2L, wbar) %*% basis$U / r_B, U_cell = contrasts(n_cells),
    T = length(calendar))
}

#' Collapsed design: one pseudo-row per month with column rows.
#'
#' With omega_it = (1 / d_it^2) / S_t, the month enters through ybar_t =
#' sum_i omega_it y_it, its cell means Bbar and Gbar and noise 1 / (kappa_T S_t).
column_pseudo_rows <- function(design) {
  rows <- design$rows
  inv_d2 <- 1 / rows$d^2
  months <- sort(unique(rows$t))
  S <- drop(rowsum(inv_d2, rows$t))
  omega <- inv_d2 / S[match(rows$t, months)]
  Bbar <- rowsum(omega * design$B_col[rows$cell, , drop = FALSE], rows$t)
  Gbar <- rowsum(omega * design$U_cell[rows$cell, , drop = FALSE], rows$t)
  dimnames(Bbar) <- dimnames(Gbar) <- NULL
  list(t = as.integer(months), S = unname(S), ybar = unname(drop(rowsum(omega * rows$y, rows$t))),
    Bbar = Bbar, Gbar = Gbar)
}

# Stan data of the column block (rows in month order, per-month start and count).
.column_stan_data <- function(design) {
  rows <- design$rows
  count <- tabulate(rows$t, nbins = design$T)
  data <- list(n_cells = design$n_cells, K_col = design$K, B_col = design$B_col,
    U_cell = design$U_cell, n_col = design$n_col, y_col = rows$y, d_col = rows$d,
    cell_col = as.integer(rows$cell), t_col = as.integer(rows$t),
    start_col = as.integer(pmax(1L, cumsum(count) - count + 1L)), count_col = as.integer(count))
  for (name in c("y_col", "d_col", "cell_col", "t_col", "start_col", "count_col"))
    dim(data[[name]]) <- length(data[[name]])
  data
}

# Kernel weights normalized per point, with the largest log kernel subtracted first.
.column_support <- function(points, knots, width_m) {
  log_kernel <- -(outer(points[, 1], knots[, 1], "-")^2 + outer(points[, 2], knots[, 2], "-")^2) /
    width_m^2
  kernel <- exp(log_kernel - apply(log_kernel, 1L, max))
  unname(kernel / rowSums(kernel))
}

# Quadrature points of each cell (bilinear map of its SW, SE, NE, NW corners) and
# their weights, which sum to one per cell.
.column_quadrature <- function(column_cells) {
  cell_ids <- unique(column_cells$cell_id)
  corner_order <- c("SW", "SE", "NE", "NW")
  index <- match(paste(rep(cell_ids, each = 4L), rep(corner_order, length(cell_ids))),
    paste(column_cells$cell_id, column_cells$corner))
  if (anyNA(index)) stop("Every column cell needs its four corners.", call. = FALSE)
  corners <- column_cells[index, , drop = FALSE]
  order <- length(gauss_legendre_3$nodes)
  n_per_cell <- order * order
  nodes <- (gauss_legendre_3$nodes + 1) / 2
  u <- rep(nodes, times = order)
  v <- rep(nodes, each = order)
  product_weight <- as.vector(outer(gauss_legendre_3$weights / 2, gauss_legendre_3$weights / 2))
  points <- vector("list", length(cell_ids))
  weight <- vector("list", length(cell_ids))
  area <- numeric(length(cell_ids))
  for (i in seq_along(cell_ids)) {
    xy <- as.matrix(corners[4L * (i - 1L) + seq_len(4L), c("x_m", "y_m")])
    origin <- xy[1, ]
    east <- xy[2, ] - origin
    north <- xy[4, ] - origin
    interaction <- (xy[3, ] - xy[4, ]) - east
    p <- outer(u, east) + outer(v, north) + outer(u * v, interaction)
    du <- matrix(east, n_per_cell, 2L, byrow = TRUE) + outer(v, interaction)
    dv <- matrix(north, n_per_cell, 2L, byrow = TRUE) + outer(u, interaction)
    jacobian <- du[, 1] * dv[, 2] - du[, 2] * dv[, 1]
    if (any(!is.finite(jacobian)) || any(jacobian <= 0))
      stop("Column cell ", cell_ids[i], " has an invalid Jacobian.", call. = FALSE)
    raw <- product_weight * jacobian
    area[i] <- sum(raw)
    weight[[i]] <- raw / area[i]
    # Translate after the geometry to keep precision in UTM metres.
    points[[i]] <- sweep(p, 2L, origin, "+")
  }
  list(points = do.call(rbind, points), weight = unlist(weight),
    cell = rep(seq_along(cell_ids), each = n_per_cell), cell_ids = cell_ids, area_m2 = area)
}
