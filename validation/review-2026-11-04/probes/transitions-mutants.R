# Mutation check, in memory only: override one internal function in the loaded
# namespace, run the area's test files, restore. Nothing on disk changes.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
Sys.setenv(NOT_CRAN = "true")
ns <- asNamespace("latents")
tdir <- "/Users/mohammedsaqr/Documents/Github/latents/tests/testthat"
run <- function(label, name, mutant, files) {
  original <- get(name, ns)
  unlockBinding(name, ns); assign(name, mutant(original), ns); lockBinding(name, ns)
  on.exit({ unlockBinding(name, ns); assign(name, original, ns); lockBinding(name, ns) })
  t0 <- proc.time()[["elapsed"]]
  res <- do.call(rbind, lapply(files, function(f) {
    r <- as.data.frame(suppressWarnings(suppressMessages(testthat::test_file(file.path(tdir, f), reporter = "silent", stop_on_failure = FALSE, env = new.env(parent = ns)))))
    data.frame(file = f, failed = sum(r$failed), errors = sum(r$error), passed = sum(r$passed))
  }))
  cat(sprintf("%-52s failed %3d  errors %2d  passed %4d  (%0.0fs) %s\n", label, sum(res$failed), sum(res$errors), sum(res$passed), proc.time()[["elapsed"]] - t0,
              if (sum(res$failed) + sum(res$errors) == 0L) "<-- SURVIVED" else ""))
}
wrap_body <- function(transform) function(f) { force(f); function(...) transform(f(...)) }
args <- commandArgs(TRUE)
sets <- list(
  a = function() {
    run("natural Jacobian: transition block x 1.1", ".multilpa_transition_jacobian",
        wrap_body(function(J) { n <- nrow(J); J[(n - 1L):n, ] <- J[(n - 1L):n, ] * 1.1; J }),
        c("test-transition-inference.R", "test-transitions-numerics.R"))
    run("prevalence: uniform over classes", ".multilpa_transition_prevalence",
        wrap_body(function(P) { P[] <- P[, ncol(P):1L]; P }),
        c("test-transitions.R", "test-release-transitions.R"))
  },
  b = function() {
    run("get_tna: class-mean instead of counts", ".multilpa_aggregate_transitions",
        function(f) { force(f); function(x) { out <- f(x); n <- x$n_profiles
          out$probabilities[] <- matrix(matrix(x$transition_probabilities, n * n, x$n_group_classes) %*% x$group_probabilities, n, n); out } },
        c("test-tna.R", "test-tna-fits.R"))
    run("chart Jacobian: identity (raw SEs wrong)", ".lta_chart_jacobian",
        function(f) { force(f); function(fit, chart, n_theta, direction = "to_raw") diag(n_theta) },
        c("test-transition-inference.R", "test-review-1026-lta-r2.R", "test-lta-general.R"))
  },
  c = function() {
    run("relabel: initial logits not rebased", ".latents_rebase_logits",
        function(f) { force(f); function(eta, order) eta[, order, drop = FALSE] },
        c("test-lta-bootstrap.R", "test-review-1026-lta.R"))
    run("second-order table: previous and from swapped", ".lta_second_order_table",
        wrap_body(function(t) { hold <- t$previous; t$previous <- t$from; t$from <- hold; t }),
        c("test-lta-general.R", "test-review-lta.R"))
  },
  d = function() {
    run("transition table: estimated flag always TRUE", ".multilpa_transition_frame",
        wrap_body(function(t) { t$estimated <- TRUE; t }),
        c("test-transitions.R", "test-review-1015-transitions.R", "test-release-transitions.R"))
    run("expected counts: multiplied by 1.01 in the table", ".multilpa_transition_frame",
        wrap_body(function(t) { t$expected_count <- t$expected_count * 1.01; t }),
        c("test-transitions.R", "test-review-1015-transitions.R", "test-release-transitions.R"))
  },
  e = function() {
    run("second-order count: K^2 replaced by K", ".lta_count_parameters",
        function(f) { force(f); function(n_profiles, n_types, d, n_categories, variance_model, n_measurement, p0, p, order = 1L, n_stayers = 0L, covariance_model = "diagonal", structure = NULL, p_second = p) {
          f(n_profiles, n_types, d, n_categories, variance_model, n_measurement, p0, p, order, n_stayers, covariance_model, structure, p_second) - (if (order >= 2L) (n_types - n_stayers) * p_second * (n_profiles - 1L) * (n_profiles^2 - n_profiles) else 0L) } },
        c("test-lta-general.R", "test-review-lta.R"))
    run("stayer: initial logits not counted", ".lta_count_parameters",
        function(f) { force(f); function(n_profiles, n_types, d, n_categories, variance_model, n_measurement, p0, p, order = 1L, n_stayers = 0L, covariance_model = "diagonal", structure = NULL, p_second = p) {
          as.integer(f(n_profiles, n_types, d, n_categories, variance_model, n_measurement, p0, p, order, n_stayers, covariance_model, structure, p_second) - n_stayers * p0 * (n_profiles - 1L)) } },
        c("test-lta-general.R", "test-review-lta.R"))
  })
wide <- list(
  f = function() run("natural Jacobian: transition block x 1.1 [wide]", ".multilpa_transition_jacobian", wrap_body(function(J) { n <- nrow(J); J[(n - 1L):n, ] <- J[(n - 1L):n, ] * 1.1; J }), c("test-transition-inference.R", "test-transitions.R")),
  g = function() run("prevalence: columns reversed [wide]", ".multilpa_transition_prevalence", wrap_body(function(P) { P[] <- P[, ncol(P):1L]; P }), c("test-transitions.R", "test-transitions-numerics.R", "test-review-1015-transitions.R", "test-lta-bootstrap.R")),
  h = function() run("chart Jacobian: identity [wide]", ".lta_chart_jacobian", function(f) { force(f); function(fit, chart, n_theta, direction = "to_raw") diag(n_theta) }, c("test-review-1026-lta-open.R", "test-audit-inference-repairs.R", "test-review-1101-bootstrap-tables.R")),
  i = function() run("estimated flag always TRUE [wide]", ".multilpa_transition_frame", wrap_body(function(t) { t$estimated <- TRUE; t }), c("test-release-transitions.R", "test-review-1015-transitions.R", "test-transitions.R", "test-tna.R", "test-review-1026-lta-r2.R")),
  j = function() run("expected counts x 1.01 [wide]", ".multilpa_transition_frame", wrap_body(function(t) { t$expected_count <- t$expected_count * 1.01; t }), c("test-tna.R", "test-transitions-numerics.R", "test-transitions.R")),
  k = function() run("stayer initial logits not counted [wide]", ".lta_count_parameters", function(f) { force(f); function(n_profiles, n_types, d, n_categories, variance_model, n_measurement, p0, p, order = 1L, n_stayers = 0L, covariance_model = "diagonal", structure = NULL, p_second = p) { as.integer(f(n_profiles, n_types, d, n_categories, variance_model, n_measurement, p0, p, order, n_stayers, covariance_model, structure, p_second) - n_stayers * p0 * (n_profiles - 1L)) } }, c("test-audit-1101-transitions-three-step.R", "test-audit-inference-repairs.R", "test-compare-lrt.R", "test-release-units.R", "test-review-1015-transitions.R", "test-review-1026-lta-r2.R", "test-review-lta.R", "test-search.R")))
if (args[1] %in% names(wide)) { wide[[args[1]]]() } else sets[[args[1]]]()
