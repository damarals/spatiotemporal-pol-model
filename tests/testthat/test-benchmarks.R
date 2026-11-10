# Empirical benchmarks (class and pooled), pool definitions, and location
# bootstrap uncertainty.

test_benchmarks <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) cached <<- score_benchmarks(test_inputs(), variants = TRUE)
    cached
  }
})

testthat::test_that("empirical_crps matches the exact EDF definition", {
  pool <- c(6, 10, 14, 18)
  obs <- c(8, 12)
  exact <- vapply(obs, function(y) mean(abs(pool - y)) -
    sum(abs(outer(pool, pool, "-"))) / (2 * length(pool)^2), numeric(1))
  testthat::expect_equal(empirical_crps(pool, obs), exact, tolerance = 1e-12)
  testthat::expect_error(empirical_crps(numeric(), 1), "non-empty")
  testthat::expect_error(empirical_crps(1, NA_real_), "finite")
})

testthat::test_that("benchmark pools and CRPS match the paper's values on every unit", {
  skip_without_moqa()
  scores <- test_benchmarks()
  b <- scores$scores
  pool_sizes <- function(u, pol, kind) {
    unique(b$n_training[b$unit == u & b$pollutant == pol & b$benchmark == kind])
  }
  # New traffic fold 1: avenue (5) and bus terminal (7); pooled without background (12 NO2, 14 PM2.5);
  # pooled with background (44 NO2, 38 PM2.5).
  testthat::expect_setequal(pool_sizes("new-traffic-1", "NO2", "class"), c(5L, 7L))
  testthat::expect_identical(pool_sizes("new-traffic-1", "NO2", "pooled"), 12L)
  testthat::expect_identical(pool_sizes("new-traffic-1", "PM25", "pooled"), 14L)
  testthat::expect_identical(pool_sizes("new-traffic-1", "NO2", "pooled_with_background"), 44L)
  testthat::expect_identical(pool_sizes("new-traffic-1", "PM25", "pooled_with_background"), 38L)
  # Omitted months at the background location: own history (26), all sites (43), without background (17).
  testthat::expect_identical(pool_sizes("omitted-months", "NO2", "class"), 26L)
  testthat::expect_identical(pool_sizes("omitted-months", "NO2", "pooled"), 43L)
  testthat::expect_identical(pool_sizes("omitted-months", "NO2", "pooled_without_background"), 17L)
  # Whole year 2023: one background PM2.5 observation before 2023 (point mass).
  testthat::expect_identical(pool_sizes("year-2023", "PM25", "class"), 1L)
  # New devices fold 1: all 13 training low-cost observations for both benchmarks.
  testthat::expect_identical(pool_sizes("new-devices-1", "PM25", "class"), 13L)
  testthat::expect_identical(pool_sizes("new-devices-1", "PM25", "pooled"), 13L)
  # Withheld background location has no class benchmark (56 targets), and new-devices benchmarks coincide.
  testthat::expect_identical(unique(paste(scores$unsupported$unit, scores$unsupported$benchmark)),
    "withheld-background class")
  testthat::expect_identical(nrow(scores$unsupported), 56L)
  dev <- b[startsWith(b$unit, "new-devices"), , drop = FALSE]
  testthat::expect_identical(dev$crps[dev$benchmark == "class"], dev$crps[dev$benchmark == "pooled"])

  # Equal-location benchmark CRPS values reported in Table 4 and Section 4.1.
  loc_crps <- function(grp, pol, kind) {
    rows <- .group_rows(b[b$benchmark == kind, , drop = FALSE], grp, pol)
    mean(tapply(rows$crps, rows$physical_location_id, mean))
  }
  testthat::expect_lt(abs(loc_crps("new-traffic", "NO2", "class") - 2.099379), 5e-7)
  testthat::expect_lt(abs(loc_crps("new-traffic", "NO2", "pooled") - 1.843377), 5e-7)
  testthat::expect_lt(abs(loc_crps("new-traffic", "PM25", "class") - 1.319997), 5e-7)
  testthat::expect_lt(abs(loc_crps("new-traffic", "PM25", "pooled") - 1.185904), 5e-7)
  testthat::expect_lt(abs(loc_crps("new-traffic", "NO2", "pooled_with_background") - 5.303086), 5e-7)
  testthat::expect_lt(abs(loc_crps("new-traffic", "PM25", "pooled_with_background") - 1.115778), 5e-7)
  testthat::expect_lt(abs(loc_crps("new-devices", "PM25", "class") - 1.892099), 5e-7)
  testthat::expect_equal(round(loc_crps("withheld-background", "NO2", "pooled"), 3), 8.570)
  testthat::expect_equal(round(loc_crps("withheld-background", "PM25", "pooled"), 3), 1.325)
  testthat::expect_equal(round(loc_crps("omitted-months", "NO2", "class"), 3), 1.676)
  testthat::expect_equal(round(loc_crps("omitted-months", "NO2", "pooled"), 3), 2.318)
  testthat::expect_equal(round(loc_crps("whole-years", "NO2", "class"), 3), 1.599)
  testthat::expect_equal(round(loc_crps("whole-years", "NO2", "pooled"), 3), 3.503)
  testthat::expect_equal(round(loc_crps("whole-years", "PM25", "class"), 3), 1.697)
  testthat::expect_equal(round(loc_crps("whole-years", "PM25", "pooled"), 3), 1.242)
})

testthat::test_that("summarise_differences computes between-location SE, stability and bootstrap intervals", {
  # Multi-location case: 3 locations, all positive and all drop-1 means positive -> stable.
  d <- c(0.4, 0.6, 0.2, 0.3, 0.5)
  loc <- c("l1", "l1", "l2", "l3", "l3")
  u <- summarise_differences(d, loc, seed = 42L, B = 500L)
  dl <- c(0.5, 0.2, 0.4)
  testthat::expect_equal(u$effect, mean(dl))
  testthat::expect_equal(u$se_between, stats::sd(dl) / sqrt(3))
  testthat::expect_identical(u$se_basis, "locations")
  testthat::expect_identical(u$agree, 3L)
  testthat::expect_true(u$stable)
  testthat::expect_lte(u$boot_lo, u$effect)
  testthat::expect_gte(u$boot_hi, u$effect)
  # One-location case: SE over targets, NA for location-level stability and bootstrap.
  u1 <- summarise_differences(d, rep("l1", 5L), seed = 42L)
  testthat::expect_equal(u1$effect, mean(d))
  testthat::expect_equal(u1$se_between, stats::sd(d) / sqrt(5))
  testthat::expect_identical(u1$se_basis, "targets")
  testthat::expect_true(is.na(u1$stable))
})