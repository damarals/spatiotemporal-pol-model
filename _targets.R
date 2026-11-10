# Pipeline from data/ to every file in results/: inputs, Stan data, the model
# fits, held-out predictions and scores, result tables, the numbers quoted in
# the text, concentration maps, predictive and TROPOMI column checks,
# simulation-based calibration and Figures 1-4. Run from the repository root
# with targets::tar_make().
# CmdStan files go under work/; the store keeps file paths and small tables.
library(targets)
tar_source(c("R", "data-raw/R/write.R"))

# Validation fits: each model with the units it is fitted on (models$units in
# R/models.R; tests/testthat/test-pipeline.R checks that the two agree).
validation_fits <- utils::read.table(header = TRUE, stringsAsFactors = FALSE, text = "
  model                                unit
  full                                 omitted-months
  full                                 new-traffic-1
  full                                 new-traffic-2
  full                                 new-traffic-3
  full                                 withheld-background
  full                                 year-2023
  full                                 year-2024
  full                                 year-2025
  full                                 new-devices-1
  full                                 new-devices-2
  full                                 new-devices-3
  full-unsigned                        new-traffic-1
  full-unsigned                        new-traffic-2
  full-unsigned                        new-traffic-3
  full-no-tropomi                      omitted-months
  full-no-tropomi                      new-traffic-1
  full-no-tropomi                      new-traffic-2
  full-no-tropomi                      new-traffic-3
  full-no-tropomi                      withheld-background
  full-no-tropomi                      year-2023
  full-no-tropomi                      year-2024
  full-no-tropomi                      year-2025
  full-no-tropomi                      new-devices-1
  full-no-tropomi                      new-devices-2
  full-no-tropomi                      new-devices-3
  full-no-tropomi-without-r            omitted-months
  full-no-tropomi-without-r            new-traffic-1
  full-no-tropomi-without-r            new-traffic-2
  full-no-tropomi-without-r            new-traffic-3
  full-no-tropomi-without-r            year-2023
  full-no-tropomi-without-r            year-2024
  full-no-tropomi-without-r            year-2025
  full-no-tropomi-without-r            new-devices-1
  full-no-tropomi-without-r            new-devices-2
  full-no-tropomi-without-r            new-devices-3
  full-no-tropomi-without-no2-profile  omitted-months
  full-no-tropomi-without-no2-profile  new-traffic-1
  full-no-tropomi-without-no2-profile  new-traffic-2
  full-no-tropomi-without-no2-profile  new-traffic-3
  full-no-tropomi-without-no2-profile  withheld-background
  full-no-tropomi-without-no2-profile  year-2023
  full-no-tropomi-without-no2-profile  year-2024
  full-no-tropomi-without-no2-profile  year-2025
  full-no-tropomi-without-low-cost     omitted-months
  full-no-tropomi-without-low-cost     new-traffic-1
  full-no-tropomi-without-low-cost     new-traffic-2
  full-no-tropomi-without-low-cost     new-traffic-3
  full-no-tropomi-without-low-cost     year-2025
  full-no-tropomi-without-pm25-profile omitted-months
  full-no-tropomi-without-pm25-profile new-traffic-1
  full-no-tropomi-without-pm25-profile new-traffic-2
  full-no-tropomi-without-pm25-profile new-traffic-3
  full-no-tropomi-without-pm25-profile withheld-background
  full-no-tropomi-without-pm25-profile year-2023
  full-no-tropomi-without-pm25-profile year-2024
  full-no-tropomi-without-pm25-profile year-2025
  full-no-tropomi-without-pm25-profile new-devices-1
  full-no-tropomi-without-pm25-profile new-devices-2
  full-no-tropomi-without-pm25-profile new-devices-3
  full-no-tropomi-without-deviation    omitted-months
  full-no-tropomi-without-deviation    new-traffic-1
  full-no-tropomi-without-deviation    new-traffic-2
  full-no-tropomi-without-deviation    new-traffic-3
  full-no-tropomi-without-deviation    withheld-background
  full-no-tropomi-without-deviation    year-2023
  full-no-tropomi-without-deviation    year-2024
  full-no-tropomi-without-deviation    year-2025
  full-no-tropomi-without-deviation    new-devices-1
  full-no-tropomi-without-deviation    new-devices-2
  full-no-tropomi-without-deviation    new-devices-3
  full-no-tropomi-without-era5         omitted-months
  full-no-tropomi-without-era5         new-traffic-1
  full-no-tropomi-without-era5         new-traffic-2
  full-no-tropomi-without-era5         new-traffic-3
  full-no-tropomi-without-era5         withheld-background
  full-no-tropomi-without-era5         year-2023
  full-no-tropomi-without-era5         year-2024
  full-no-tropomi-without-era5         year-2025
  full-no-tropomi-free-gain            omitted-months
  full-no-tropomi-free-gain            new-traffic-1
  full-no-tropomi-free-gain            new-traffic-2
  full-no-tropomi-free-gain            new-traffic-3
  full-no-tropomi-free-gain            withheld-background
  full-no-tropomi-free-gain            year-2023
  full-no-tropomi-free-gain            year-2024
  full-no-tropomi-free-gain            year-2025
  full-no-tropomi-without-seasonal     omitted-months
  full-no-tropomi-without-seasonal     new-traffic-1
  full-no-tropomi-without-seasonal     new-traffic-2
  full-no-tropomi-without-seasonal     new-traffic-3
  full-no-tropomi-without-seasonal     withheld-background
  full-no-tropomi-without-seasonal     year-2023
  full-no-tropomi-without-seasonal     year-2024
  full-no-tropomi-without-seasonal     year-2025
  full-no-tropomi-without-satpm25      omitted-months
  full-no-tropomi-without-satpm25      new-traffic-1
  full-no-tropomi-without-satpm25      new-traffic-2
  full-no-tropomi-without-satpm25      new-traffic-3
  full-no-tropomi-without-satpm25      withheld-background
  full-no-tropomi-without-satpm25      year-2023
  full-no-tropomi-without-satpm25      year-2024
  full-no-tropomi-without-satpm25      year-2025
  full-no-tropomi-without-amplitude    omitted-months
  full-no-tropomi-without-amplitude    new-traffic-1
  full-no-tropomi-without-amplitude    new-traffic-2
  full-no-tropomi-without-amplitude    new-traffic-3
  full-no-tropomi-without-amplitude    withheld-background
  full-no-tropomi-without-amplitude    year-2023
  full-no-tropomi-without-amplitude    year-2024
  full-no-tropomi-without-amplitude    year-2025
  full-no-tropomi-one-class            omitted-months
  full-no-tropomi-one-class            new-traffic-1
  full-no-tropomi-one-class            new-traffic-2
  full-no-tropomi-one-class            new-traffic-3
  full-no-tropomi-one-class            withheld-background
  full-no-tropomi-one-class            year-2023
  full-no-tropomi-one-class            year-2024
  full-no-tropomi-one-class            year-2025
  full-no-tropomi-without-profile-cor  new-traffic-1
  full-no-tropomi-without-profile-cor  new-traffic-2
  full-no-tropomi-without-profile-cor  new-traffic-3
  full-no-tropomi-without-profile-cor  withheld-background
  full-no-tropomi-without-no2-cov      omitted-months
  full-no-tropomi-without-no2-cov      new-traffic-1
  full-no-tropomi-without-no2-cov      new-traffic-2
  full-no-tropomi-without-no2-cov      new-traffic-3
  full-no-tropomi-without-no2-cov      withheld-background
  full-no-tropomi-without-no2-cov      year-2023
  full-no-tropomi-without-no2-cov      year-2024
  full-no-tropomi-without-no2-cov      year-2025
  full-no-tropomi-without-pm25-cov     new-traffic-1
  full-no-tropomi-without-pm25-cov     new-traffic-2
  full-no-tropomi-without-pm25-cov     new-traffic-3
  full-no-tropomi-without-pm25-cov     withheld-background
  full-no-tropomi-without-pm25-cov     new-devices-1
  full-no-tropomi-without-pm25-cov     new-devices-2
  full-no-tropomi-without-pm25-cov     new-devices-3
  reduced                              omitted-months
  reduced                              new-traffic-1
  reduced                              new-traffic-2
  reduced                              new-traffic-3
  reduced                              withheld-background
  reduced                              year-2023
  reduced                              year-2024
  reduced                              year-2025
  reduced                              new-devices-1
  reduced                              new-devices-2
  reduced                              new-devices-3
  reduced-free-gain                    omitted-months
  reduced-free-gain                    new-traffic-1
  reduced-free-gain                    new-traffic-2
  reduced-free-gain                    new-traffic-3
  reduced-free-gain                    withheld-background
  reduced-free-gain                    year-2023
  reduced-free-gain                    year-2024
  reduced-free-gain                    year-2025
  reduced-without-no2-cov              omitted-months
  reduced-without-no2-cov              new-traffic-1
  reduced-without-no2-cov              new-traffic-2
  reduced-without-no2-cov              new-traffic-3
  reduced-without-no2-cov              withheld-background
  reduced-without-no2-cov              year-2023
  reduced-without-no2-cov              year-2024
  reduced-without-no2-cov              year-2025
  reduced-without-pm25-profile         omitted-months
  reduced-without-pm25-profile         new-traffic-1
  reduced-without-pm25-profile         new-traffic-2
  reduced-without-pm25-profile         new-traffic-3
  reduced-without-pm25-profile         withheld-background
  reduced-without-pm25-profile         year-2023
  reduced-without-pm25-profile         year-2024
  reduced-without-pm25-profile         year-2025
  reduced-without-pm25-profile         new-devices-1
  reduced-without-pm25-profile         new-devices-2
  reduced-without-pm25-profile         new-devices-3
  reduced-without-low-cost             omitted-months
  reduced-without-low-cost             new-traffic-1
  reduced-without-low-cost             new-traffic-2
  reduced-without-low-cost             new-traffic-3
  reduced-without-low-cost             year-2025
  reduced-without-deviation            omitted-months
  reduced-without-deviation            new-traffic-1
  reduced-without-deviation            new-traffic-2
  reduced-without-deviation            new-traffic-3
  reduced-without-deviation            withheld-background
  reduced-without-deviation            year-2023
  reduced-without-deviation            year-2024
  reduced-without-deviation            year-2025
  reduced-without-deviation            new-devices-1
  reduced-without-deviation            new-devices-2
  reduced-without-deviation            new-devices-3
  reduced-one-class                    omitted-months
  reduced-one-class                    new-traffic-1
  reduced-one-class                    new-traffic-2
  reduced-one-class                    new-traffic-3
  reduced-one-class                    withheld-background
  reduced-one-class                    year-2023
  reduced-one-class                    year-2024
  reduced-one-class                    year-2025
  final-no-tropomi                     omitted-months
  final-no-tropomi                     new-traffic-1
  final-no-tropomi                     new-traffic-2
  final-no-tropomi                     new-traffic-3
  final-no-tropomi                     withheld-background
  final-no-tropomi                     year-2023
  final-no-tropomi                     year-2024
  final-no-tropomi                     year-2025
  final-no-tropomi                     new-devices-1
  final-no-tropomi                     new-devices-2
  final-no-tropomi                     new-devices-3
  final-no-tropomi-without-low-cost    omitted-months
  final-no-tropomi-without-low-cost    new-traffic-1
  final-no-tropomi-without-low-cost    new-traffic-2
  final-no-tropomi-without-low-cost    new-traffic-3
  final-no-tropomi-without-low-cost    year-2025
  final                                omitted-months
  final                                new-traffic-1
  final                                new-traffic-2
  final                                new-traffic-3
  final                                withheld-background
  final                                year-2023
  final                                year-2024
  final                                year-2025
  final                                new-devices-1
  final                                new-devices-2
  final                                new-devices-3
  final-one-class                      omitted-months
  final-one-class                      new-traffic-1
  final-one-class                      new-traffic-2
  final-one-class                      new-traffic-3
  final-eight-knots                    omitted-months
  final-eight-knots                    new-traffic-1
  final-eight-knots                    new-traffic-2
  final-eight-knots                    new-traffic-3
  final-eight-knots                    withheld-background
  final-half-normal                    omitted-months
  final-half-normal                    new-traffic-1
  final-half-normal                    new-traffic-2
  final-half-normal                    new-traffic-3
  final-half-normal                    withheld-background
")

list(
  # Inputs (only the files the pipeline reads, data_files() in R/inputs.R) and the Stan program.
  tar_target(data_paths, data_files("data"), format = "file"),
  tar_target(input, read_inputs(files = data_paths)),
  tar_target(grid, read_grid(files = data_paths)),
  tar_target(program_file, "stan/fusion.stan", format = "file"),
  tar_target(program, compile_program("work/stan", program_file), format = "file"),

  # Fits: one branch per model and unit, then the final model and the final model
  # without TROPOMI on all data 2020-2025.
  tar_target(fit_pairs, validation_fits),
  tar_target(fit, run_fit(input, fit_pairs$unit, fit_pairs$model, program),
    pattern = map(fit_pairs), format = "file"),
  tar_target(final_fit_files, run_fit(input, "all-data", "final", program), format = "file"),
  tar_target(final_fit, read_fit(input, "all-data", "final", dirname(final_fit_files[1L]))),
  tar_target(no_tropomi_fit_files, run_fit(input, "all-data", "final-no-tropomi", program),
    format = "file"),
  tar_target(no_tropomi_fit, read_fit(input, "all-data", "final-no-tropomi",
    dirname(no_tropomi_fit_files[1L]))),

  # Held-out scores, benchmarks and sampler diagnostics.
  tar_target(scores, score_fit(input, fit_pairs$unit, fit_pairs$model, dirname(fit[1L])),
    pattern = map(fit_pairs, fit)),
  tar_target(fit_diagnostics,
    diagnostics_table(list(read_fit(input, fit_pairs$unit, fit_pairs$model, dirname(fit[1L])))),
    pattern = map(fit_pairs, fit)),
  tar_target(benchmarks, score_benchmarks(input, variants = TRUE)),

  # Tables.
  tar_target(scores_csv, write_result(scores_table(scores, benchmarks), "results/scores.csv"),
    format = "file"),
  tar_target(elimination_csv, write_result(elimination_table(scores), "results/elimination.csv"),
    format = "file"),
  tar_target(full_vs_final_csv, write_result(full_vs_final_table(scores),
    "results/full-vs-final.csv"), format = "file"),
  tar_target(benchmarks_csv, write_result(compare_benchmarks(benchmarks, scores,
    models = c("final", "final-no-tropomi", "full-no-tropomi", "reduced", "full"), variants = TRUE),
    "results/benchmarks.csv"), format = "file"),
  tar_target(posterior_csv, write_result(posterior_table(final_fit), "results/posterior.csv"),
    format = "file"),
  tar_target(diagnostics_csv, write_result(diagnostics_table(rbind(fit_diagnostics,
    diagnostics_table(list(final_fit, no_tropomi_fit)))), "results/diagnostics.csv"), format = "file"),

  # Numbers quoted in the text: per-target scores, coverage, the background
  # location by year, gains by location, the pooled TROPOMI contrast, derived
  # posterior summaries and the calibration factor of the low-cost network.
  tar_target(targets_csv, write_result(target_table(scores), "results/targets.csv"), format = "file"),
  tar_target(coverage_csv, write_result(coverage_table(scores), "results/coverage.csv"),
    format = "file"),
  tar_target(background_csv, write_result(background_table(scores, benchmarks),
    "results/background.csv"), format = "file"),
  tar_target(location_gains_csv, write_result(location_gains_table(scores),
    "results/location-gains.csv"), format = "file"),
  tar_target(pooled_contrast_csv, write_result(pooled_contrast(scores[scores$model == "final", ],
    scores[scores$model == "final-no-tropomi", ]), "results/pooled-contrast.csv"), format = "file"),
  tar_target(derived_csv, write_result(derived_table(list(final_fit, no_tropomi_fit)),
    "results/derived.csv"), format = "file"),
  tar_target(calibration_csv, write_result(calibration_factor(final_fit, input),
    "results/calibration-factor.csv"), format = "file"),

  # Maps and checks of the fits to all data.
  tar_target(maps, compute_maps(final_fit, input, grid)),
  tar_target(map_table, map_cells(maps, grid)),
  tar_target(map_cells_csv, write_table(map_table, "results/map-cells.csv"), format = "file"),
  tar_target(no_tropomi_maps, compute_maps(no_tropomi_fit, input, grid)),
  tar_target(map_summary_csv, write_result(rbind(map_summary(maps, grid, "final"),
    map_summary(no_tropomi_maps, grid, "final-no-tropomi")), "results/map-summary.csv"),
    format = "file"),
  tar_target(ppc_csv, write_table(public_ppc(posterior_predictive_checks(final_fit)), "results/ppc.csv"),
    format = "file"),
  tar_target(column_check, column_checks(final_fit, input)),
  tar_target(column_checks_csv, write_result(column_check$summary, "results/column-checks.csv"),
    format = "file"),

  # Simulation-based calibration: replication i uses seed sbc_seed + i.
  tar_target(sbc_design, prepare_fit(input, partition(input, "all-data"), "final")),
  tar_target(sbc_replications, seq_len(120L)),
  tar_target(sbc_replication, sbc_replicate(sbc_replications, sbc_design,
    file.path("work", "sbc", sbc_replications), program),
    pattern = map(sbc_replications), iteration = "list"),
  tar_target(sbc_reps, sbc_table(sbc_replication)),
  tar_target(sbc_summary, summarise_sbc(sbc_reps)),
  tar_target(sbc_csv, write_table(sbc_summary, "results/sbc.csv"), format = "file"),
  tar_target(sbc_replications_csv, write_table(sbc_reps, "results/sbc-replications.csv"),
    format = "file"),
  tar_target(sbc_bins_csv, write_result(sbc_bins(sbc_reps), "results/sbc-bins.csv"), format = "file"),
  tar_target(sbc_overview_csv, write_result(sbc_overview(sbc_reps, sbc_summary),
    "results/sbc-overview.csv"), format = "file"),

  # Figures 1-4.
  tar_target(figure_1, plot_study_area("data", "results/figures/study-area.pdf", data_paths),
    format = "file"),
  tar_target(figure_2, plot_source_coverage(input, "results/figures/source-coverage.pdf"),
    format = "file"),
  tar_target(figure_3, plot_column_checks(column_check$summary, "results/figures/column-checks.pdf"),
    format = "file"),
  tar_target(figure_4, plot_maps(map_table, input$sites, data_file(data_paths, "boundary"),
    "results/figures/maps.pdf"), format = "file")
)
