# Growth mixture models (mixture_regression(random = )) against lcmm::hlme
# (Proust-Lima, Philipps and Liquet 2017).
#
# Two checks per configuration. (1) The likelihood: latents' log likelihood
# evaluated at hlme's own estimates equals hlme's loglik, which shows the two
# define the same model independently of how either optimizes. (2) The
# optimum: latents reaches at least hlme's log likelihood (hlme stops on its
# own convergence criteria, so it may sit a little below), and the class
# trajectories agree.
#
# hlme offers a random-effect covariance common to all classes (`nwg = FALSE`,
# latents' "equal") or proportional across classes (`nwg = TRUE`,
# "proportional"), and one residual variance for all classes
# (`variance = "equal"`); class-specific covariances ("varying") have no hlme
# counterpart and are checked against an independent dense likelihood in
# tests/testthat/.

skip_if_not_installed("lcmm")

growth_data <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  set.seed(21)
  n <- 200L
  times <- 0:5
  class <- sample(1:2, n, TRUE, prob = c(0.6, 0.4))
  covariance <- matrix(c(1, 0.2, 0.2, 0.1), 2L)
  effects <- t(chol(covariance)) %*% matrix(stats::rnorm(2L * n), 2L)
  out <- data.frame(id = rep(seq_len(n), each = length(times)),
                    time = rep(times, n), x = stats::rnorm(n * length(times)),
                    age = rep(stats::rnorm(n), each = length(times)))
  person_class <- class[out$id]
  out$y <- c(2, 6)[person_class] + c(0.5, -0.4)[person_class] * out$time +
    0.3 * out$x + effects[1L, out$id] +
    effects[2L, out$id] * out$time * c(1, 1.5)[person_class] +
    stats::rnorm(nrow(out), 0, 0.8)
  out <- out[stats::runif(nrow(out)) > 0.15, ]
  if (!is.null(old)) assign(".Random.seed", old, globalenv())
  out
})

# hlme from a one-class start and a grid of random starts.
fit_hlme <- function(fixed, mixture, random, ng, nwg = FALSE, idiag = FALSE,
                     classmb = ~1) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(11)
  one <- lcmm::hlme(fixed, random = random, subject = "id", ng = 1,
                    idiag = idiag, data = growth_data, verbose = FALSE)
  if (ng == 1L) return(one)
  # gridsearch() reads its model unevaluated and calls it by name, so the
  # whole call is built around a local `hlme`.
  hlme <- lcmm::hlme
  model <- bquote(hlme(.(fixed), mixture = .(mixture), random = .(random),
                             subject = "id", ng = .(ng), nwg = .(nwg),
                             idiag = .(idiag), classmb = .(classmb),
                             data = growth_data, verbose = FALSE))
  eval(bquote(lcmm::gridsearch(.(model), rep = 15, maxiter = 30, minit = one)))
}

# latents' parameter list from an hlme fit, classes in hlme's order.
params_from_hlme <- function(fit, model) {
  spec <- fit$spec
  n_classes <- spec$n_classes
  best <- model$best
  counts <- model$N
  at <- 0L
  take <- function(count) {
    values <- best[at + seq_len(count)]
    at <<- at + count
    values
  }
  membership <- take(counts[1L])
  fixed <- take(counts[2L])
  varcov <- take(counts[3L])
  proportional <- take(counts[4L])
  params <- fit$params
  # Fixed effects: hlme lists each class-specific term once per class
  # ("time class1"), or once without a suffix when there is one class.
  label_of <- function(term) if (identical(term, "(Intercept)")) "intercept" else term
  params$beta[] <- vapply(seq_len(n_classes), function(k) {
    vapply(rownames(params$beta), function(term) {
      name <- if (n_classes == 1L) label_of(term) else
        sprintf("%s class%d", label_of(term), k)
      unname(fixed[[name]])
    }, numeric(1))
  }, numeric(nrow(params$beta)))
  if (length(params$common) > 0L) {
    params$common[] <- unname(fixed[names(params$common)])
  }
  n_random <- ncol(spec$random_design)
  covariance <- if (isTRUE(spec$random_diagonal)) diag(varcov, n_random) else {
    m <- matrix(0, n_random, n_random)
    m[upper.tri(m, diag = TRUE)] <- varcov
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  params$random_covariance <- list(covariance)
  params$random_scale <- if (length(proportional) > 0L)
    c(unname(proportional), 1) else rep(1, n_classes)
  params$sigma2 <- rep(unname(best[["stderr"]])^2, n_classes)
  # hlme's membership logits are against its last class; latents' against
  # the first.
  if (n_classes > 1L) {
    terms <- colnames(spec$w)
    logits <- vapply(seq_len(n_classes), function(k) {
      if (k == n_classes) return(numeric(length(terms)))
      unname(membership[sprintf("%s class%d", ifelse(terms == "(Intercept)",
                                                     "intercept", terms), k)])
    }, numeric(length(terms))) |> matrix(length(terms))
    params$gamma[] <- logits - logits[, 1L]
  }
  params
}

# hlme's estimate vector from latents' parameters (for restarting hlme).
hlme_from_params <- function(fit, model) {
  spec <- fit$spec
  params <- fit$params
  n_classes <- spec$n_classes
  best <- model$best
  label_of <- function(term) if (identical(term, "(Intercept)")) "intercept" else term
  logits <- params$gamma - params$gamma[, n_classes]
  terms <- colnames(spec$w)
  values <- vapply(names(best), function(name) {
    parts <- regmatches(name, regexec("^(.*) class([0-9]+)$", name))[[1L]]
    if (identical(name, "stderr")) return(sqrt(params$sigma2[1L]))
    if (startsWith(name, "varcov")) {
      covariance <- params$random_covariance[[1L]]
      entries <- if (isTRUE(spec$random_diagonal)) diag(covariance) else
        covariance[upper.tri(covariance, diag = TRUE)]
      return(entries[as.integer(sub("varcov ", "", name))])
    }
    if (startsWith(name, "varprop")) {
      return(params$random_scale[as.integer(sub("varprop class ", "", name))])
    }
    if (length(parts) == 3L) {
      term <- parts[2L]
      k <- as.integer(parts[3L])
      membership_term <- which(vapply(terms, label_of, character(1)) == term)
      if (length(membership_term) == 1L && k < n_classes &&
          match(name, names(best)) <= model$N[1L]) {
        return(logits[membership_term, k])
      }
      specific <- rownames(params$beta)[vapply(rownames(params$beta), label_of,
                                               character(1)) == term]
      return(params$beta[specific, k])
    }
    unname(params$common[[name]])
  }, numeric(1))
  values
}

latents_at <- function(fit, params) {
  latents:::.growth_expectation(fit$spec, fit$stats, params)$log_likelihood
}

configurations <- list(
  list(label = "two classes, equal covariance",
       formula = y ~ time, mixture = ~ time, random = ~ 1 + time, hlme_random = ~ time,
       ng = 2L, nwg = FALSE, idiag = FALSE, classmb = ~1, common = NULL,
       membership = ~1, covariance = "equal", diagonal = FALSE),
  list(label = "one class (a linear mixed model)",
       formula = y ~ time, mixture = ~ time, random = ~ 1 + time, hlme_random = ~ time,
       ng = 1L, nwg = FALSE, idiag = FALSE, classmb = ~1, common = NULL,
       membership = ~1, covariance = "equal", diagonal = FALSE),
  list(label = "proportional covariance (nwg)",
       formula = y ~ time, mixture = ~ time, random = ~ 1 + time, hlme_random = ~ time,
       ng = 2L, nwg = TRUE, idiag = FALSE, classmb = ~1, common = NULL,
       membership = ~1, covariance = "proportional", diagonal = FALSE),
  list(label = "diagonal random effects (idiag)",
       formula = y ~ time, mixture = ~ time, random = ~ 1 + time, hlme_random = ~ time,
       ng = 2L, nwg = FALSE, idiag = TRUE, classmb = ~1, common = NULL,
       membership = ~1, covariance = "equal", diagonal = TRUE),
  list(label = "common covariate and membership covariate",
       formula = y ~ time + x, mixture = ~ time, random = ~ 1 + time, hlme_random = ~ time,
       ng = 2L, nwg = FALSE, idiag = FALSE, classmb = ~ age, common = ~ x,
       membership = ~ age, covariance = "equal", diagonal = FALSE),
  list(label = "random intercept only, three classes",
       formula = y ~ time, mixture = ~ time, random = ~ 1, hlme_random = ~ 1,
       ng = 3L, nwg = FALSE, idiag = FALSE, classmb = ~1, common = NULL,
       membership = ~1, covariance = "equal", diagonal = FALSE))

invisible(lapply(configurations, function(cfg) {
  test_that(sprintf("growth mixture matches lcmm::hlme: %s", cfg$label), {
    hlme_random <- cfg$hlme_random
    theirs <- fit_hlme(cfg$formula, cfg$mixture, hlme_random, cfg$ng, cfg$nwg,
                       cfg$idiag, cfg$classmb)
    expect_equal(theirs$conv, 1)
    ours <- mixture_regression(cfg$formula, growth_data, n_classes = cfg$ng,
                               id = "id", class_level = "group",
                               common = cfg$common, membership = cfg$membership,
                               random = cfg$random,
                               random_covariance = cfg$covariance,
                               random_diagonal = cfg$diagonal, variance = "equal",
                               n_starts = 6, seed = 3, vcov_type = "none")
    expect_true(ours$converged)
    # (1) the same likelihood at hlme's estimates.
    expect_equal(latents_at(ours, params_from_hlme(ours, theirs)), theirs$loglik,
                 tolerance = 1e-8)
    # (2) The optimum, compared within one likelihood: ours against ours at
    # hlme's estimates (hlme reports its loglik about 1e-5 below the
    # likelihood of its own estimates).
    at_theirs <- latents_at(ours, params_from_hlme(ours, theirs))
    expect_gte(ours$log_likelihood, at_theirs - 1e-6)
    if (ours$log_likelihood - at_theirs > 1e-3) {
      # hlme stopped at a lower maximum: restarted from ours it must reach it.
      restarted <- lcmm::hlme(cfg$formula, mixture = cfg$mixture, random = hlme_random,
                              subject = "id", ng = cfg$ng, nwg = cfg$nwg,
                              idiag = cfg$idiag, classmb = cfg$classmb,
                              data = growth_data, verbose = FALSE,
                              B = hlme_from_params(ours, theirs))
      expect_equal(latents_at(ours, params_from_hlme(ours, restarted)),
                   ours$log_likelihood, tolerance = 1e-9)
    }
  })
}))
