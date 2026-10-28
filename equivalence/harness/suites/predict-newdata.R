# predict() on new data, against the reference programs' own scoring code.
#
# The full-output suites compare in-sample posteriors only. Scoring rows a fit
# never saw goes through a separate code path (`.multilpa_prediction_rows()`:
# column checks, categorical label matching, missing indicators under FIML,
# new groups classified from their own rows; `.mixture_new_spec()`: the fitted
# design rebuilt on new data), so it needs its own evidence.
#
# Design. Each fixture is split deterministically: every fifth row (or every
# fifth group) is held out, the package is fitted to the rest, and the held-out
# rows are predicted with every `type` the method offers. The reference is the
# other program's scoring code evaluated at the package's own estimates. One
# formula at one fixed parameter set is compared, so nothing but the two
# implementations of that formula separates the numbers. Whether the package's
# estimates are the reference's maximum is a separate claim, made by the
# mclust, poLCA, glca and flexmix suites on the full data.
#
#   Gaussian, faithful, full covariance, varying (mclust "VVV", 2 components).
#     mclust's `estep()` gives the posteriors and `dens(logarithm = TRUE)` the
#     log mixture density. Mapping: `pro` <-> `profile_probabilities[1, ]`,
#     `mean` <-> `t(means)`, `variance$sigma` <-> `covariances`
#     (indicator x indicator x profile), with `cholsigma` its Cholesky
#     factors. The fit uses `missing = "fiml"` and one held-out row lacks
#     `waiting`; its reference is mclust's univariate "V" model at the
#     observed indicator's marginal (the mean and variance of that
#     coordinate), which is what integrating out a missing Gaussian
#     coordinate means.
#   Categorical, poLCA `values` (4 binary items, 2 classes).
#     `poLCA.posterior(lc, y)`, poLCA's own E-step in C, with `lc$probs[[j]]`
#     <-> `response_probabilities[[j]]` (profiles x categories, columns in
#     `categorical_levels` order), `lc$P` <-> `profile_probabilities[1, ]`
#     and `y` the category codes, missing coded 0 (poLCA skips it). The log
#     density is `log sum_k P_k f_k(y)` with `f_k` from poLCA's own
#     `poLCA.ylik.C()`, unscaled as `poLCA()` unscales it (it returns
#     `f_k * .Machine$double.xmax`; leaving that in was the one disagreement
#     the first run of this suite showed, a reference slip of exactly
#     log(double.xmax) = 709.78). One held-out row lacks `C` (fit under FIML).
#   Two-level categorical, the package's `student_esm` (3 binary activities,
#     2 profiles x 2 group classes, every fifth student held out). There is no
#     external program that scores new groups, so the reference is the
#     two-level likelihood written from its definition in base R
#     (`two_level_evaluate()` in two-level-lca-reference.R, which the glca
#     suite checks against glca's own likelihood and posteriors). Mapping:
#     `delta` <-> `group_probabilities`, `gamma` <-> `profile_probabilities`,
#     `rho` <-> `response_probabilities`. The row density predict() reports is
#     that of a row from a group not yet seen: one group class whose profile
#     shares are `delta %*% gamma`, each row its own group. One held-out row
#     lacks an activity (fit under FIML).
#   Regression mixture, the flexmix suite's `gaussian-2-membership` fixture
#     (y ~ x, membership ~ z, 2 classes). flexmix's own `posterior()`,
#     `prior()` and `predict()` on `newdata`, with the flexmix object's
#     components rebuilt by its driver's `defineComponent()` at the package's
#     estimates, as the flexmix suite's `.flexmix_at_estimates()` does:
#     component `coef` <-> `params$beta[, k]`, `sigma` <-> sqrt(`params$sigma2`),
#     concomitant `coef` <-> `params$gamma` (class 1 the zero baseline in
#     both). `type = "response"` is `sum_k prior_k mean_k`.
#
# Not compared, deliberately:
#   * predict() for membership-covariate fits (`multilpa_covariates`) and
#     latent transition fits (`multilpa_transitions`) is not implemented: both
#     refuse with `latents_unsupported_prediction`, as `?predict.multilpa`
#     documents. The suite pins that refusal. If either is implemented, the
#     refusal rows fail and the reference named in the review plan (EQ-10)
#     has to be written: `poLCA.posterior(lc, y, x)` for covariates, and a
#     depmixS4 model on the new sequences with `setpars()` and
#     `forwardbackward()$gamma` for transitions.
#   * `type = "probabilities"` of a regression mixture exists for the ordinal
#     family only, which flexmix cannot fit; for this Gaussian fixture its
#     refusal (`latents_bad_argument`) is checked.
#   * The fits themselves are not put in a ledger here: every output of fits of
#     these kinds is compared in the mclust, polca, glca and flexmix suites.
#
# Tolerances (relative, `|r - o| / max(|r|, 1)`, as `compare_block()` uses):
#   * 1e-10 for every posterior, prior, class mean, fitted value and log
#     density: the same formula at the same parameters, so only floating-point
#     order of operations separates the two (measured 1e-16 to 4e-15 when the
#     suite was written). Modal classes and group identifiers are integers or
#     labels and are therefore compared exactly in effect.
#   * 1e-8 for the log density summed over the held-out rows, a quantity of
#     tens to hundreds whose rounding error grows with the number of rows.
#
# Refusals of new rows the fit cannot score (an unseen category, a missing
# column, a missing predictor) are rows outside the ledger, as in the
# `refusals` suite: each passes when the expected condition class is raised.

source(file.path("equivalence", "harness", "two-level-lca-reference.R"), local = environment())

suite_predict_newdata <- function() {
  require_suite_packages(
    "mclust", "poLCA", "flexmix",
    reason = paste("the Gaussian, categorical and regression-mixture scoring code",
                   "predict() is checked against, and the faithful and values data"))
  # mclust's estep() and dens() dispatch to their model-specific functions by
  # name, so the namespace has to be attached rather than only loaded.
  if (!"package:mclust" %in% search()) {
    suppressPackageStartupMessages(library(mclust))
    on.exit(detach("package:mclust", unload = FALSE), add = TRUE, after = FALSE)
  }
  stack_cases(list(.predict_gaussian_case(), .predict_categorical_case(),
                   .predict_two_level_case(), .predict_regression_case(),
                   .predict_unsupported_case()))
}

# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------

#' The held-out rows or groups: every fifth one
#' @param n Number of rows or groups.
#' @return A logical vector, `TRUE` for the held-out ones.
#' @noRd
.predict_held_out <- function(n) seq_len(n) %% 5L == 0L

#' Run code under a seed and restore the caller's random stream
#' @param seed Integer seed.
#' @param code Code to evaluate.
#' @return The value of `code`.
#' @noRd
.predict_with_seed <- function(seed, code) {
  saved <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(saved)) {
    if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", saved, globalenv()), add = TRUE, after = FALSE)
  set.seed(seed)
  code
}

#' The tolerances of this suite (see the header)
.predict_tolerance <- c(default = 1e-10)

#' One refusal row: the call must raise the expected condition class
#'
#' The condition is recorded, not swallowed: its classes go into the row's
#' label, and anything but the expected class (another class, an unclassed
#' error, or a result) fails the row.
#'
#' @param label Description of the case.
#' @param call A function of no arguments making the call.
#' @param expected The condition class it must raise.
#' @return A comparison row.
#' @noRd
.predict_refusal <- function(label, call, expected) {
  seen <- tryCatch({
    call()
    "a result"
  }, error = function(condition) paste(class(condition), collapse = "/"))
  passed <- expected %in% strsplit(seen, "/", fixed = TRUE)[[1L]]
  compare_values(sprintf("%s: refused with %s (got: %s)", label, expected, seen),
                 reference = 1, obtained = as.numeric(passed), tolerance = 0.5)
}

#' The three tables predict() gives for a mixture of profiles or classes
#' @param fit A `multilpa` fit.
#' @param newdata The held-out rows.
#' @return A named list of keyed tables.
#' @noRd
.predict_profile_tables <- function(fit, newdata) {
  list(predict_class = keyed_table(predict(fit, newdata, type = "class"), "row"),
       predict_posterior = keyed_table(predict(fit, newdata, type = "posterior"), "row"),
       predict_density = keyed_table(predict(fit, newdata, type = "density"), "row"))
}

#' Expected predict() tables from reference posteriors and log densities
#'
#' The modal profile takes the first on a tie, as `?predict.multilpa`
#' documents. A two-level fit adds the group identifier and the modal group
#' class with its posterior.
#'
#' @param posteriors Rows-by-profiles reference posteriors.
#' @param log_density Each row's reference log mixture density.
#' @param groups `NULL`, or a list with `id` (the identifier column name),
#'   `values` (each row's group) and `posteriors` (rows-by-group-classes, each
#'   row's group's posterior).
#' @return A named list of tables shaped like [.predict_profile_tables()]'s.
#' @noRd
.predict_expected_profiles <- function(posteriors, log_density, groups = NULL) {
  n <- nrow(posteriors)
  modal <- max.col(posteriors, ties.method = "first")
  base <- data.frame(row = seq_len(n))
  if (!is.null(groups)) base[[groups$id]] <- groups$values
  base$profile <- modal
  base$posterior <- posteriors[cbind(seq_len(n), modal)]
  if (!is.null(groups)) {
    group_modal <- max.col(groups$posteriors, ties.method = "first")
    base$group_class <- group_modal
    base$group_posterior <- groups$posteriors[cbind(seq_len(n), group_modal)]
  }
  probabilities <- stats::setNames(as.data.frame(posteriors),
                                   sprintf("probability_profile_%d", seq_len(ncol(posteriors))))
  list(predict_class = base,
       predict_posterior = cbind(base, probabilities),
       predict_density = data.frame(row = seq_len(n), log_density = log_density))
}

#' Compare the predict() tables and the summed log density, and label the rows
#' @param label Case label.
#' @param obtained,expected Named lists of tables.
#' @param source,dataset Report labels.
#' @return A `harness_case`.
#' @noRd
.predict_compare <- function(label, obtained, expected, source, dataset, extra = NULL) {
  summed <- compare_values(
    sprintf("%s: log density summed over the %d held-out rows", label,
            nrow(expected$predict_density)),
    sum(expected$predict_density$log_density), sum(obtained$predict_density$log_density),
    tolerance = 1e-8, scale = "relative")
  case <- compare_verb(label, obtained, expected, tolerance = .predict_tolerance,
                       extra = rbind(.with_output_columns(summed),
                                     if (!is.null(extra)) .with_output_columns(extra)))
  case$comparisons <- transform(case$comparisons, source = source, dataset = dataset,
                                precision = "machine")
  case
}

#' Label refusal rows for the report
#' @param rows Comparison rows from [.predict_refusal()].
#' @param dataset Report label.
#' @return The rows with `source`, `dataset` and `precision`.
#' @noRd
.predict_refusal_rows <- function(rows, dataset) {
  transform(.with_output_columns(rows), source = "refusal parity", dataset = dataset,
            precision = "machine")
}

# ---------------------------------------------------------------------------
# Gaussian: mclust
# ---------------------------------------------------------------------------

#' mclust's posteriors and log densities at given parameters
#'
#' Rows are scored by missingness pattern: the observed coordinates' marginal
#' is the Gaussian with those coordinates' means and covariance block, which
#' mclust scores as "VVV" (two or more coordinates) or "V" (one).
#'
#' @param x Rows-by-indicators matrix, `NA` where missing.
#' @param shares Mixing proportions.
#' @param means Indicators-by-components matrix.
#' @param sigma Indicators x indicators x components array.
#' @return A list of `posteriors` and `log_density`.
#' @noRd
.predict_mclust_score <- function(x, shares, means, sigma) {
  patterns <- apply(!is.na(x), 1L, paste, collapse = "")
  scored <- lapply(unique(patterns), function(pattern) {
    rows <- which(patterns == pattern)
    observed <- !is.na(x[rows[[1L]], ])
    data <- x[rows, observed, drop = FALSE]
    g <- length(shares)
    parameters <- if (sum(observed) == 1L) {
      list(pro = shares, mean = means[observed, ],
           variance = list(modelName = "V", d = 1L, G = g,
                           sigmasq = sigma[observed, observed, ]))
    } else {
      block <- sigma[observed, observed, , drop = FALSE]
      list(pro = shares, mean = means[observed, , drop = FALSE],
           variance = list(modelName = "VVV", d = sum(observed), G = g, sigma = block,
                           cholsigma = array(vapply(seq_len(g), function(k) chol(block[, , k]),
                                                    block[, , 1L]), dim(block))))
    }
    model <- parameters$variance$modelName
    data <- if (identical(model, "V")) drop(data) else data
    list(rows = rows,
         posteriors = matrix(mclust::estep(data = data, modelName = model,
                                           parameters = parameters)$z, length(rows)),
         log_density = mclust::dens(data = data, modelName = model,
                                    parameters = parameters, logarithm = TRUE))
  })
  rows <- unlist(lapply(scored, `[[`, "rows"))
  posteriors <- do.call(rbind, lapply(scored, `[[`, "posteriors"))
  log_density <- unlist(lapply(scored, `[[`, "log_density"))
  stopifnot("every row must be scored once" = identical(sort(rows), seq_len(nrow(x))))
  list(posteriors = posteriors[order(rows), , drop = FALSE],
       log_density = log_density[order(rows)])
}

#' faithful, full covariance with varying volumes, against mclust
#' @return A `harness_case`.
#' @noRd
.predict_gaussian_case <- function() {
  frame <- datasets::faithful
  indicators <- c("eruptions", "waiting")
  held_out <- .predict_held_out(nrow(frame))
  fit <- lpa(frame[!held_out, ], indicators, 2L, covariance_model = "full",
             variance_model = "varying", missing = "fiml", n_starts = 10L, seed = 1L,
             tol = 1e-12)
  newdata <- frame[held_out, indicators]
  row.names(newdata) <- NULL
  # A held-out row with an indicator missing, scored under FIML.
  newdata$waiting[[3L]] <- NA
  reference <- .predict_mclust_score(
    as.matrix(newdata), shares = fit$profile_probabilities[1L, ],
    means = t(fit$means[, indicators]), sigma = fit$covariances[indicators, indicators, ])
  label <- "faithful, mclust VVV, 2 components, 54 held-out rows"
  refusals <- .predict_refusal(sprintf("%s: newdata lacking `waiting`", label),
                               \() predict(fit, newdata["eruptions"]), "latents_bad_data")
  case <- .predict_compare(label, .predict_profile_tables(fit, newdata),
                           .predict_expected_profiles(reference$posteriors, reference$log_density),
                           source = sprintf("mclust %s", utils::packageVersion("mclust")),
                           dataset = "faithful")
  stack_cases(list(case, .predict_refusal_rows(refusals, "faithful")))
}

# ---------------------------------------------------------------------------
# Categorical: poLCA
# ---------------------------------------------------------------------------

#' values, two classes, against poLCA's E-step
#' @return A `harness_case`.
#' @noRd
.predict_categorical_case <- function() {
  values <- get(utils::data(list = "values", package = "poLCA", envir = environment()),
                envir = environment())
  items <- c("A", "B", "C", "D")
  held_out <- .predict_held_out(nrow(values))
  fit <- lca(values[!held_out, items], items, 2L, missing = "fiml", n_starts = 10L,
             seed = 1L, tol = 1e-12)
  newdata <- values[held_out, items]
  row.names(newdata) <- NULL
  newdata$C[[2L]] <- NA
  # poLCA's coding: categories 1..K in the fit's category order, 0 missing.
  codes <- vapply(items, function(item) {
    match(as.character(newdata[[item]]), fit$categorical_levels[[item]])
  }, integer(nrow(newdata)))
  codes[is.na(codes)] <- 0L
  # poLCA.posterior() reads `probs`, `P` and `x` (one column: no covariates).
  at_estimates <- list(probs = unname(fit$response_probabilities[items]),
                       P = fit$profile_probabilities[1L, ], x = matrix(1, 1L, 1L))
  posteriors <- poLCA::poLCA.posterior(at_estimates, codes)
  # poLCA.ylik.C() returns each class's likelihood multiplied by
  # .Machine$double.xmax (against underflow); poLCA()'s own log-likelihood
  # expression takes the scale back out, and so does this one.
  conditional <- poLCA:::poLCA.ylik.C(poLCA:::poLCA.vectorize(at_estimates$probs), codes)
  log_density <- log(drop(conditional %*% at_estimates$P)) - log(.Machine$double.xmax)
  label <- "values, poLCA, 2 classes, 43 held-out rows"
  unseen <- newdata
  unseen$A <- as.character(unseen$A)
  unseen$A[[1L]] <- "3"
  refusals <- .predict_refusal(sprintf("%s: a category the fit never saw", label),
                               \() predict(fit, unseen), "latents_bad_data")
  case <- .predict_compare(label, .predict_profile_tables(fit, newdata),
                           .predict_expected_profiles(posteriors, log_density),
                           source = sprintf("poLCA %s", utils::packageVersion("poLCA")),
                           dataset = "values")
  stack_cases(list(case, .predict_refusal_rows(refusals, "values")))
}

# ---------------------------------------------------------------------------
# Two-level categorical: the base-R two-level likelihood
# ---------------------------------------------------------------------------

#' student_esm, 2 profiles x 2 group classes, new students
#' @return A `harness_case`.
#' @noRd
.predict_two_level_case <- function() {
  activities <- c("time_with_friends", "on_social_media", "tv_video_games")
  esm <- student_esm[c("student", activities)]
  students <- unique(esm$student)
  held_students <- students[.predict_held_out(length(students))]
  held_out <- esm$student %in% held_students
  fit <- multilca(esm[!held_out, ], activities, "student", 2L, 2L, missing = "fiml",
                  n_starts = 10L, seed = 1L, tol = 1e-12)
  newdata <- esm[held_out, ]
  row.names(newdata) <- NULL
  newdata$on_social_media[[3L]] <- NA
  codes <- .two_level_codes(newdata, fit$categorical_levels[activities])
  rho <- fit$response_probabilities[activities]
  scored <- two_level_evaluate(fit$group_probabilities, fit$profile_probabilities, rho,
                               codes, newdata$student)
  # A row's density as a member of an unseen group: the group classes'
  # profile shares averaged over the group classes, each row its own group.
  shares <- matrix(drop(fit$group_probabilities %*% fit$profile_probabilities), 1L)
  alone <- two_level_evaluate(1, shares, rho, codes, seq_len(nrow(newdata)))
  expected <- .predict_expected_profiles(
    scored$posteriors, unname(alone$group_log_likelihood),
    groups = list(id = "student", values = newdata$student,
                  posteriors = scored$group_posteriors[as.character(newdata$student), ,
                                                       drop = FALSE]))
  label <- sprintf("student_esm, 2 profiles x 2 group classes, %d held-out students (%d rows)",
                   length(held_students), nrow(newdata))
  unseen <- newdata
  unseen$tv_video_games <- as.character(unseen$tv_video_games)
  unseen$tv_video_games[[1L]] <- "maybe"
  refusals <- rbind(
    .predict_refusal(sprintf("%s: a category the fit never saw", label),
                     \() predict(fit, unseen), "latents_bad_data"),
    .predict_refusal(sprintf("%s: newdata lacking the group identifier", label),
                     \() predict(fit, newdata[activities]), "latents_bad_data"))
  case <- .predict_compare(label, .predict_profile_tables(fit, newdata), expected,
                           source = "base-R two-level likelihood (two-level-lca-reference.R)",
                           dataset = "student_esm")
  stack_cases(list(case, .predict_refusal_rows(refusals, "student_esm")))
}

# ---------------------------------------------------------------------------
# Regression mixture: flexmix
# ---------------------------------------------------------------------------

#' The flexmix suite's gaussian-2-membership data, regenerated
#'
#' The same draws as `.flexmix_fixtures()` makes for that fixture (seeds 501
#' and 5011), so the held-out rows are a fifth of the data that suite fits.
#'
#' @return A data frame with `x`, `z`, `r`, `truth` and `y`.
#' @noRd
.predict_regression_data <- function() {
  frame <- .predict_with_seed(501L, {
    frame <- data.frame(x = stats::rnorm(400L), z = stats::rnorm(400L), r = stats::rnorm(400L))
    frame$truth <- 1L + (stats::runif(nrow(frame)) < stats::plogis(0.3 + 0.9 * frame$z))
    frame
  })
  .predict_with_seed(5011L, {
    frame$y <- ifelse(frame$truth == 1L, 1 + 2 * frame$x, -1 - frame$x) +
      stats::rnorm(nrow(frame), sd = ifelse(frame$truth == 1L, 0.7, 1))
  })
  frame
}

#' A flexmix object at the package's estimates
#'
#' flexmix is run for one EM step from the true classes, only to obtain an
#' object of the right model; its components are then rebuilt by the driver's
#' own `defineComponent()` from the package's coefficients and residual SDs,
#' and its concomitant coefficients replaced by the package's membership
#' logits.
#'
#' @param fit A `latents_mixture_regression` fit.
#' @param train The rows it was fitted to.
#' @return A `flexmix` object.
#' @noRd
.predict_flexmix_at <- function(fit, train) {
  scaffold <- flexmix::flexmix(y ~ x, data = train, k = 2L, cluster = train$truth,
                               concomitant = flexmix::FLXPmultinom(~ z),
                               control = list(iter.max = 1L, minprior = 0, verbose = 0))
  driver <- scaffold@model[[1L]]
  at <- scaffold
  at@components <- lapply(seq_len(2L), function(class) {
    component <- scaffold@components[[class]][[1L]]
    parameters <- component@parameters
    stopifnot("every flexmix coefficient must have a package estimate" =
                all(names(parameters$coef) %in% rownames(fit$params$beta)))
    parameters$coef[] <- fit$params$beta[names(parameters$coef), class]
    parameters$sigma <- sqrt(fit$params$sigma2[[class]])
    list(driver@defineComponent(c(parameters, list(df = component@df))))
  })
  stopifnot("the membership terms must be flexmix's concomitant terms, in order" =
              identical(rownames(at@concomitant@coef), colnames(fit$spec$w)),
            "class 1 must be the zero baseline in both" =
              all(at@concomitant@coef[, 1L] == 0) && all(fit$params$gamma[, 1L] == 0))
  at@concomitant@coef[] <- fit$params$gamma
  at
}

#' gaussian-2-membership, against flexmix's scoring of new rows
#' @return A `harness_case`.
#' @noRd
.predict_regression_case <- function() {
  frame <- .predict_regression_data()
  held_out <- .predict_held_out(nrow(frame))
  train <- frame[!held_out, ]
  fit <- mixture_regression(y ~ x, train, 2L, membership = ~ z, n_starts = 10L,
                            seed = 501L, tol = 1e-12, max_iter = 20000L)
  newdata <- frame[held_out, c("x", "z", "y")]
  row.names(newdata) <- NULL
  at <- .predict_flexmix_at(fit, train)
  prior <- unname(flexmix::prior(at, newdata = newdata))
  posteriors <- unname(flexmix::posterior(at, newdata = newdata))
  means <- unname(do.call(cbind, flexmix::predict(at, newdata = newdata)))
  n <- nrow(newdata)
  modal <- max.col(posteriors, ties.method = "first")
  expected <- list(
    predict_response = data.frame(row = seq_len(n), fitted = rowSums(prior * means)),
    predict_class_response = data.frame(row = rep(seq_len(n), 2L),
                                        class = rep(c("class_1", "class_2"), each = n),
                                        prior = as.vector(prior), fitted = as.vector(means)),
    predict_posterior = data.frame(row = seq_len(n), class = modal,
                                   posterior = posteriors[cbind(seq_len(n), modal)],
                                   probability_class_1 = posteriors[, 1L],
                                   probability_class_2 = posteriors[, 2L]))
  obtained <- list(
    predict_response = keyed_table(predict(fit, newdata, type = "response"), "row"),
    predict_class_response = keyed_table(predict(fit, newdata, type = "class_response"),
                                         c("row", "class")),
    predict_posterior = keyed_table(predict(fit, newdata, type = "posterior"), "row"))
  label <- sprintf("gaussian-2-membership, flexmix, 2 classes, %d held-out rows", n)
  missing_predictor <- newdata
  missing_predictor$x[[2L]] <- NA
  refusals <- rbind(
    .predict_refusal(sprintf("%s: type = \"probabilities\" for a Gaussian outcome", label),
                     \() predict(fit, newdata, type = "probabilities"), "latents_bad_argument"),
    .predict_refusal(sprintf("%s: a missing predictor", label),
                     \() predict(fit, missing_predictor), "latents_bad_data"))
  case <- compare_verb(label, obtained, expected, tolerance = .predict_tolerance)
  case$comparisons <- transform(case$comparisons,
                                source = sprintf("flexmix %s", utils::packageVersion("flexmix")),
                                dataset = "gaussian-2-membership", precision = "machine")
  stack_cases(list(case, .predict_refusal_rows(refusals, "gaussian-2-membership")))
}

# ---------------------------------------------------------------------------
# Families predict() does not score yet
# ---------------------------------------------------------------------------

#' Membership-covariate and transition fits refuse new data by class
#'
#' Small simulated fits, held-out rows or sequences offered to predict(). The
#' transition data follow the depmixS4 suite's generator (stay probability
#' 0.75, state means -2 and 2) at a smaller size.
#'
#' @return A comparison frame.
#' @noRd
.predict_unsupported_case <- function() {
  classes <- .predict_with_seed(20261027L, {
    n <- 200L
    frame <- data.frame(x = stats::rnorm(n))
    truth <- 1L + (stats::runif(n) < stats::plogis(0.8 * frame$x))
    items <- vapply(seq_len(3L), function(item) {
      ifelse(stats::runif(n) < ifelse(truth == 1L, 0.2, 0.8), "yes", "no")
    }, character(n))
    cbind(frame, stats::setNames(as.data.frame(items), c("u1", "u2", "u3")))
  })
  held_out <- .predict_held_out(nrow(classes))
  covariate_fit <- lca(classes[!held_out, ], c("u1", "u2", "u3"), 2L,
                       profile_covariates = "x", n_starts = 2L, seed = 1L)
  sequences <- .predict_with_seed(20261028L, {
    groups <- 30L
    occasions <- 5L
    transition <- matrix(c(0.75, 0.25, 0.25, 0.75), 2L)
    frame <- do.call(rbind, lapply(seq_len(groups), function(group) {
      path <- unlist(Reduce(function(previous, step) {
        sample.int(2L, 1L, prob = transition[previous, ])
      }, seq_len(occasions - 1L), init = sample.int(2L, 1L), accumulate = TRUE))
      data.frame(g = group, t = seq_len(occasions), state = path)
    }))
    frame$y1 <- stats::rnorm(nrow(frame), c(-2, 2)[frame$state], 0.8)
    frame$y2 <- stats::rnorm(nrow(frame), c(-1.5, 1.5)[frame$state], 0.9)
    frame
  })
  held_sequences <- sequences$g %in% which(.predict_held_out(30L))
  transition_fit <- lta(sequences[!held_sequences, ], c("y1", "y2"), "g", n_profiles = 2L,
                        time = "t", n_starts = 2L, seed = 1L)
  rows <- rbind(
    .predict_refusal("membership-covariate fit (lca, profile_covariates = \"x\"), 40 held-out rows",
                     \() predict(covariate_fit, classes[held_out, ]),
                     "latents_unsupported_prediction"),
    .predict_refusal("latent transition fit (lta, 2 states), 6 held-out sequences",
                     \() predict(transition_fit, sequences[held_sequences, ]),
                     "latents_unsupported_prediction"))
  transform(.with_output_columns(rows), source = "refusal parity",
            dataset = "simulated", precision = "machine")
}
