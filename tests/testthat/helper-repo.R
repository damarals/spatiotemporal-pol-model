# Paths, inputs and the compiled program shared by the tests (working directory tests/testthat).

repo_root <- normalizePath(file.path("..", ".."))

skip_without_moqa <- function() {
  testthat::skip_if_not(file.exists(file.path(repo_root, "data", "moqa", "moqa-monthly.csv")),
    "the low-cost network data are not in data/moqa/")
}

skip_without_cmdstan <- function() {
  testthat::skip_if(!nzchar(Sys.getenv("CMDSTAN")), "CMDSTAN is not set")
}

test_inputs <- local({
  inputs <- NULL
  function() {
    if (is.null(inputs)) inputs <<- read_inputs(file.path(repo_root, "data"))
    inputs
  }
})

test_program <- local({
  exe <- NULL
  function() {
    if (is.null(exe)) exe <<- compile_program(file.path(repo_root, "work", "stan"),
      file.path(repo_root, "stan", "fusion.stan"))
    exe
  }
})

# One model on one unit, prepared as fit_model() prepares it.
test_case <- function(unit, model) {
  input <- test_inputs()
  part <- partition(input, unit)
  list(input = input, part = part, prep = prepare_fit(input, part, model), model = model_spec(model))
}
