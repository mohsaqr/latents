# Builds tests/equivalence/fixtures/lta-structure-references.rds: the retained
# datasets and the independent reference standard errors for the ten
# constrained covariance structures in the homogeneous lta()
# (`multilpa_transitions`). Run from the project root:
#
#   Rscript equivalence/lta-structure-se/make-references.R
#
# The reference is tests/equivalence/helper-lta-structure-reference.R: its own
# forward recursion and Gaussian density, the structure charts of
# helper-structure-reference.R, logits against the first category, the
# central second difference of helper-structure-reference.R for the Hessian
# (h = 2e-4) and numDeriv's Richardson Jacobian for the delta method. latents
# supplies only the point at which the curvature is taken (its estimates),
# and the reference's own gradient there is stored so that the test can
# require the point to be a maximum of the reference likelihood as well.
# Fixed before the first build: standard errors must agree within 1e-5
# relative, the log likelihood within 1e-8 relative.

devtools::load_all(".", quiet = TRUE)
source(file.path("tests", "equivalence", "helper-structure-reference.R"))
source(file.path("tests", "equivalence", "helper-lta-structure-reference.R"))
stopifnot(requireNamespace("numDeriv", quietly = TRUE))

panel <- function(n, waves, stay, seed, d = 2L, missing_share = 0, ragged = FALSE) {
  set.seed(seed)
  type <- if (length(stay) > 1L) sample.int(length(stay), n, replace = TRUE) else rep(1L, n)
  frame <- do.call(rbind, lapply(seq_len(n), function(person) {
    path <- Reduce(function(previous, wave) {
      if (stats::runif(1L) < stay[type[person]]) previous else 3L - previous
    }, seq_len(waves - 1L), init = sample.int(2L, 1L), accumulate = TRUE)
    data.frame(person = person, wave = seq_len(waves), state = path)
  }))
  s <- frame$state
  frame$a <- stats::rnorm(nrow(frame), c(0, 2)[s], c(1, 1.4)[s])
  frame$b <- stats::rnorm(nrow(frame), c(0, 1.8)[s], c(0.8, 1.1)[s]) + 0.3 * frame$a
  if (d >= 3L) frame$c <- stats::rnorm(nrow(frame), c(1, -0.5)[s], c(0.9, 0.7)[s]) -
    0.25 * frame$b
  if (missing_share > 0) frame$a[stats::runif(nrow(frame)) < missing_share] <- NA
  if (ragged) frame <- frame[!(frame$wave == waves & frame$person %% 4L == 0L), ]
  frame$state <- NULL
  rownames(frame) <- NULL
  frame
}

constrained <- c("EII", "VII", "EVI", "VEI", "EEV", "VEV", "EVE", "VVE", "VEE", "EVV")
cases <- list(
  one_class = list(data = panel(150L, 4L, 0.8, 101L), vars = c("a", "b"),
                   n_group_classes = 1L, missing = "error", structures = constrained),
  two_classes = list(data = panel(300L, 6L, c(0.95, 0.5), 102L), vars = c("a", "b"),
                     n_group_classes = 2L, missing = "error",
                     structures = c("VEI", "EVE", "VEE")),
  fiml_ragged = list(data = panel(150L, 4L, 0.8, 103L, missing_share = 0.1, ragged = TRUE),
                     vars = c("a", "b"), n_group_classes = 1L, missing = "fiml",
                     structures = c("EVI", "VVE", "EVV")),
  three_indicators = list(data = panel(160L, 4L, 0.8, 104L, d = 3L),
                          vars = c("a", "b", "c"), n_group_classes = 1L,
                          missing = "error", structures = c("VEV", "EVE", "VEE")))

fit_case <- function(case, structure) {
  suppressWarnings(lta(case$data, case$vars, "person", 2L, time = "wave",
                       n_group_classes = case$n_group_classes, model = structure,
                       missing = case$missing, n_starts = 4L, seed = 1L,
                       tol = 1e-12, max_iter = 20000L))
}

references <- lapply(names(cases), function(name) {
  case <- cases[[name]]
  results <- lapply(case$structures, function(structure) {
    fit <- fit_case(case, structure)
    stopifnot(fit$converged, !isTRUE(fit$boundary))
    reference <- reference_lta_structure_se(
      fit, case$data, jacobian = function(f, x) numDeriv::jacobian(f, x))
    chart <- reference_lta_chart(fit)
    ordered <- order(case$data$person, case$data$wave)
    sequences <- unname(split(ordered, case$data$person[ordered]))
    gradient <- numDeriv::grad(function(psi) {
      p <- chart$decode(psi)
      reference_lta_loglik(as.matrix(case$data[case$vars]), sequences, p$means,
                           p$sigma, p$initial, p$transition, p$group_probabilities)
    }, chart$start)
    message(sprintf("%s %s: log likelihood %.6f, %d coordinates", name, structure,
                    reference$log_likelihood, length(chart$start)))
    list(structure = structure, n_parameters = fit$n_parameters,
         log_likelihood = reference$log_likelihood, estimate = reference$estimate,
         standard_error = reference$standard_error, gradient = gradient,
         latents_estimate_point = list(
           means = unname(fit$means), variances = unname(fit$variances),
           covariances = if (is.null(fit$covariances)) NULL else unname(fit$covariances),
           initial = unname(fit$initial_probabilities),
           transition = unname(fit$transition_probabilities),
           group = unname(fit$group_probabilities)))
  })
  names(results) <- case$structures
  list(data = case$data, vars = case$vars, n_group_classes = case$n_group_classes,
       missing = case$missing, results = results)
})
names(references) <- names(cases)

fixture <- list(
  source = paste(
    "equivalence/lta-structure-se/make-references.R:",
    "tests/equivalence/helper-lta-structure-reference.R (independent forward",
    "recursion, Helmert/Givens/matrix-log charts, first-category logits),",
    "central second-difference Hessian (h = 2e-4), numDeriv Richardson Jacobian"),
  built = format(Sys.Date()), r_version = R.version.string,
  numderiv_version = as.character(utils::packageVersion("numDeriv")),
  latents_version = as.character(utils::packageVersion("latents")),
  fit_settings = "lta(n_profiles = 2, n_starts = 4, seed = 1, tol = 1e-12, max_iter = 20000)",
  tolerance = c(standard_error = 1e-5, log_likelihood = 1e-8),
  cases = references)
saveRDS(fixture, file.path("tests", "equivalence", "fixtures",
                           "lta-structure-references.rds"), version = 2)
