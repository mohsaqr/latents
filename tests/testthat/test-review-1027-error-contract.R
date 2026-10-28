# Review 2026-10-27: the error contract of the exported verbs. Every common
# mistake is refused with a `latents_*` class, a condition that carries no
# call (so never the data), and a message that names the argument or column
# at fault and, when there is one, the remedy. Before this review 23 of these
# 56 cases raised a bare `simpleError`, and a misspelt column, a missing
# column, a repeated column and `id` among `vars` shared one sentence.

contract_data <- function() {
  withr::local_seed(20261027)
  data.frame(id = rep(seq_len(60L), each = 4L), t = rep(1:4, 60L),
             y = stats::rnorm(240L), z = stats::rnorm(240L),
             a = factor(sample(3L, 240L, TRUE)), b = factor(sample(3L, 240L, TRUE)),
             x = stats::rnorm(240L))
}

contract_cases <- function() {
  d <- contract_data()
  fit <- suppressWarnings(lpa(d, c("y", "z"), n_profiles = 2, n_starts = 1, seed = 1))
  grid <- suppressWarnings(enumerate_lpa(d, c("y", "z"), n_profiles = 1:2, n_starts = 1,
                                         seed = 1))
  mr <- suppressWarnings(mixture_regression(y ~ x, d, n_classes = 2, n_starts = 1,
                                            seed = 1))
  env <- environment()
  case <- function(expr, class, pattern) {
    list(expr = substitute(expr), class = class, pattern = pattern, env = env)
  }
  list(
    lpa_data_fit = case(lpa(fit, c("y", "z"), n_profiles = 2), "latents_bad_data",
                        "`fit` is a fitted model.*Pass the data frame"),
    lpa_data_grid = case(lpa(grid, c("y", "z"), n_profiles = 2), "latents_bad_data",
                         "`grid` is an enumeration result"),
    multilca_data_grid = case(multilca(grid, c("a", "b"), id = "id", n_profiles = 2),
                              "latents_bad_data", "`grid` is an enumeration result"),
    lta_data_fit = case(lta(fit, c("y", "z"), id = "id", n_profiles = 2, time = "t"),
                        "latents_bad_data", "`fit` is a fitted model"),
    mixreg_data_fit = case(mixture_regression(y ~ x, fit, n_classes = 2),
                           "latents_bad_data", "`fit` is a fitted model"),
    enumerate_data_fit = case(enumerate_lpa(fit, c("y", "z"), n_profiles = 1:2),
                              "latents_bad_data", "`fit` is a fitted model"),
    vars_numeric = case(lpa(d, vars = 1, n_profiles = 2), "latents_bad_argument", "`vars`"),
    vars_typo = case(lpa(d, vars = c("y", "Z"), n_profiles = 2), "latents_bad_data",
                     "`Z`.*Did you mean `z`"),
    vars_absent = case(lpa(d, vars = c("y", "nope"), n_profiles = 2), "latents_bad_data",
                       "`nope`"),
    vars_repeated = case(lpa(d, vars = c("y", "y"), n_profiles = 2), "latents_bad_data",
                         "`y` more than once"),
    factor_in_lpa = case(lpa(d, c("y", "a"), n_profiles = 2), "latents_bad_data",
                         "`a` is not numeric.*lca\\(\\)"),
    missing_values = case(lpa(transform(d, y = replace(y, 1:3, NA)), c("y", "z"),
                              n_profiles = 2), "latents_bad_data",
                          "`y` has missing values \\(y: 3\\).*missing = \"fiml\""),
    infinite_values = case(lpa(transform(d, y = replace(y, 1, Inf)), c("y", "z"),
                               n_profiles = 2), "latents_bad_data",
                           "`y` has infinite or NaN values \\(y: 1\\)"),
    categorical_missing = case(lca(transform(d, a = replace(a, 1:2, NA)), c("a", "b"),
                                   n_classes = 2), "latents_bad_data",
                               "`a` has missing values \\(a: 2\\)"),
    id_absent = case(multilpa(d, c("y", "z"), id = "ID", n_profiles = 2), "latents_bad_data",
                     "`ID`.*Did you mean `id`"),
    id_in_vars = case(multilpa(d, c("y", "z"), id = "y", n_profiles = 2), "latents_bad_data",
                      "also listed in `vars`"),
    id_two = case(multilpa(d, c("y", "z"), id = c("id", "t"), n_profiles = 2),
                  "latents_bad_argument", "`id`"),
    time_absent = case(multilpa(d, c("y", "z"), id = "id", n_profiles = 2, time = "T"),
                       "latents_bad_data", "`T`.*Did you mean `t`"),
    lta_time_absent = case(lta(d, c("y", "z"), id = "id", n_profiles = 2, time = "T"),
                           "latents_bad_data", "`T`.*Did you mean `t`"),
    formula_text = case(mixture_regression("y ~ x", d, n_classes = 2), "latents_bad_argument",
                        "`formula`"),
    enumerate_k_text = case(enumerate_lpa(d, c("y", "z"), n_profiles = "1:2"),
                            "latents_bad_argument", "`n_profiles`"),
    three_step_data = case(three_step(d, d, "x"), "latents_not_a_fit",
                           "data frame with 240 rows"),
    three_step_outcome = case(three_step(fit, d, "nope"), "latents_bad_data", "`nope`"),
    r3step_data = case(r3step(d, d, "x"), "latents_not_a_fit", "data frame"),
    sensitivity_data = case(sensitivity(d, seeds = 1:2), "latents_not_a_fit", "`d`"),
    sensitivity_grid = case(sensitivity(grid, seeds = 1:2), "latents_not_a_fit",
                            "get_candidate\\(\\)"),
    get_candidate_fit = case(get_candidate(fit, 2), "latents_not_a_fit",
                             "already a single model"),
    candidate_fit_data = case(candidate_fit(d, 2), "latents_not_a_fit",
                              "Enumerate models"),
    diagnostics_data = case(diagnostics(d), "latents_not_a_fit", "`d`"),
    report_data = case(report(d), "latents_not_a_fit", "`d`"),
    inference_data = case(parameter_inference(d), "latents_not_a_fit", "`d`"),
    bootstrap_data = case(bootstrap_lrt(d, fit), "latents_not_a_fit", "`d`"),
    get_tna_data = case(get_tna(d), "latents_not_a_fit", "`d`"),
    get_group_tna_data = case(get_group_tna(d), "latents_not_a_fit", "`d`"),
    descriptives_by = case(descriptives(fit, "y"), "latents_bad_argument",
                           "not \"y\".*descriptives\\(data"),
    descriptives_absent = case(descriptives(d, "nope"), "latents_bad_data", "`nope`"),
    plot_diagnostics_data = case(plot_diagnostics(d), "latents_not_a_fit",
                                 "`d` is a data frame.*plot_diagnostics\\(fit\\)"),
    plot_profiles_grid = case(plot_profiles(grid), "latents_not_a_fit",
                              "`grid` is an enumeration.*plot_enumeration\\(\\)"),
    # Release review (REG-2, REG-6, REG-7, REG-8, API-1): the regression surface.
    mixreg_confint_level = case(confint(mr, level = 2), "latents_bad_argument", "`level`"),
    mixreg_results_level = case(get_results(mr, "classes", level = c(0.9, 0.95)),
                                "latents_bad_argument", "`level`"),
    mixreg_simulate_nsim = case(simulate(mr, nsim = 0), "latents_bad_argument", "`nsim`"),
    mixreg_recovery_truth = case(get_results(mr, "recovery", data = d, truth = "nope"),
                                 "latents_bad_data", "`truth`"),
    mixreg_common_intercept = case(mixture_regression(y ~ x, d, n_classes = 2,
                                                      common = "intercept"),
                                   "latents_bad_argument", "`common` names no predictor"),
    mixreg_single_level = case(mixture_regression(y ~ x + a, transform(d, w = 1 * (a == "1")),
                                                  n_classes = 1, weights = "w"),
                               "latents_bad_data", "`a` takes a single value"),
    mixreg_all_zero_weights = case(mixture_regression(y ~ x + a, transform(d, w = 0),
                                                      n_classes = 1, weights = "w"),
                                   "latents_bad_weights", "`w` must be positive"),
    enumerate_regressions_k = case(enumerate_regressions(y ~ x, d, n_classes = 0),
                                   "latents_bad_argument", "`n_classes`"),
    enumerate_regressions_bootstrap = case(enumerate_regressions(y ~ x, d, n_classes = 1:2,
                                                                 bootstrap = -1),
                                           "latents_bad_argument", "`bootstrap`"),
    enumerate_regressions_unknown = case(enumerate_regressions(y ~ x, d, n_classes = 1:2,
                                                               foo = 1),
                                         "latents_bad_argument", "does not use `foo`"))
}

test_that("every common mistake is refused with a class, no call, and a naming message", {
  cases <- contract_cases()
  invisible(lapply(names(cases), function(name) {
    case <- cases[[name]]
    condition <- tryCatch(suppressWarnings(eval(case$expr, case$env)), error = identity)
    expect_true(inherits(condition, "error"), label = name)
    expect_true(inherits(condition, case$class),
                label = sprintf("%s raises %s (got %s)", name, case$class,
                                paste(class(condition), collapse = "/")))
    expect_null(conditionCall(condition), label = name)
    expect_match(conditionMessage(condition), case$pattern, label = name)
  }))
})

test_that("a refusal of the wrong object is also a bad argument, for older handlers", {
  data <- contract_data()
  condition <- expect_error(sensitivity(data, seeds = 1:2), class = "latents_not_a_fit")
  expect_s3_class(condition, "latents_bad_argument")
})

test_that("the classed stopifnot keeps its contract", {
  expect_null(.latents_require("ok" = TRUE, "also" = c(TRUE, TRUE)))
  expect_error(.latents_require("first" = TRUE, "second fails" = FALSE),
               "second fails", class = "latents_bad_argument")
  expect_error(.latents_require("NA is not TRUE" = NA), class = "latents_bad_argument")
  expect_error(.latents_require("empty is not TRUE" = logical()),
               class = "latents_bad_argument")
  expect_error(.latents_require(1 > 2), "`1 > 2` is not TRUE")
  # Sequential: a later check is not evaluated once an earlier one fails.
  expect_error(.latents_require("stops here" = FALSE, "never" = stop("evaluated")),
               "stops here", class = "latents_bad_argument")
  expect_error(.latents_require("custom" = FALSE, .class = "latents_bad_data"),
               class = "latents_bad_data")
})

test_that("the message helpers name and suggest", {
  expect_identical(.latents_name_list("a"), "`a`")
  expect_identical(.latents_name_list(c("a", "b", "c")), "`a`, `b` and `c`")
  expect_identical(.latents_did_you_mean("Score", c("score", "id")), " Did you mean `score`?")
  expect_identical(.latents_did_you_mean("scroe", c("score", "id")), " Did you mean `score`?")
  expect_identical(.latents_did_you_mean("zzzzzz", c("score", "id")), "")
  expect_match(.latents_describe(data.frame(a = 1:3)), "3 rows and 1 column")
  expect_match(.latents_describe(1:4), "integer vector of length 4")
  expect_identical(.latents_describe(NULL), "NULL")
})

test_that("missing values are refused with per-indicator counts, and FIML only when it helps", {
  values <- data.frame(y = c(1, NA, NA), z = c(NA, 2, 3), w = 1:3)
  condition <- expect_error(.latents_refuse_missing_values(values), class = "latents_bad_data")
  expect_match(conditionMessage(condition), "`y` and `z` have missing values \\(y: 2, z: 1\\)")
  expect_match(conditionMessage(condition), "missing = \"fiml\"")
  expect_null(.latents_refuse_missing_values(values, missing = "fiml"))
  condition <- expect_error(.latents_refuse_missing_values(values, fiml_available = FALSE),
                            class = "latents_bad_data")
  expect_no_match(conditionMessage(condition), "fiml")
  # Infinite and NaN values are refused even under FIML, which treats only NA as missing.
  expect_error(.latents_refuse_missing_values(data.frame(y = c(1, Inf)), missing = "fiml"),
               "infinite or NaN", class = "latents_bad_data")
})
