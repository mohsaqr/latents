# Full-covariance transition simulation and the transition bootstrap.

boot_lta_data <- function(n = 120L, occasions = 4L, seed = 41L) {
  set.seed(seed)
  transition <- matrix(c(0.85, 0.15, 0.15, 0.85), 2L)
  states <- matrix(NA_integer_, n, occasions)
  states[, 1L] <- sample.int(2L, n, replace = TRUE)
  invisible(lapply(seq_len(occasions)[-1L], function(t) {
    states[, t] <<- vapply(states[, t - 1L], function(from)
      sample.int(2L, 1L, prob = transition[from, ]), integer(1))
  }))
  s <- as.vector(t(states))
  y1 <- c(0, 3)[s] + stats::rnorm(length(s))
  data.frame(id = rep(seq_len(n), each = occasions),
             time = rep(seq_len(occasions), n),
             y1 = y1, y2 = c(0.5, 3.5)[s] + 0.6 * y1 + stats::rnorm(length(s)))
}

boot_lta_fit <- function(data, ...) {
  suppressWarnings(lta(data, c("y1", "y2"), "id", 2, time = "time", n_starts = 3,
                       seed = 1, ...))
}

test_that("diagonal draws keep the package's random stream", {
  skip_on_cran()
  block <- list(means = matrix(c(0, 3, 1, -1), 2L, dimnames = list(NULL, c("a", "b"))),
                variances = matrix(c(1, 4, 0.5, 2), 2L, dimnames = list(NULL, c("a", "b"))))
  profile <- c(1L, 2L, 2L, 1L, 2L)
  set.seed(11)
  drawn <- .lta_draw_continuous(block, profile, c("a", "b"))
  set.seed(11)
  first <- stats::rnorm(5L, block$means[profile, "a"], sqrt(block$variances[profile, "a"]))
  second <- stats::rnorm(5L, block$means[profile, "b"], sqrt(block$variances[profile, "b"]))
  expect_identical(unname(drawn), unname(cbind(first, second)))
  expect_identical(colnames(drawn), c("a", "b"))
})

test_that("full-covariance draws reproduce each profile's mean and covariance", {
  skip_on_cran()
  covariance <- array(c(1, 0.6, 0.2, 0.6, 2, -0.5, 0.2, -0.5, 1.5,
                        0.5, -0.3, 0, -0.3, 1, 0.4, 0, 0.4, 3), c(3L, 3L, 2L))
  means <- matrix(c(0, 5, 1, -2, 3, 0), 2L, dimnames = list(NULL, c("a", "b", "c")))
  block <- list(means = means, variances = t(apply(covariance, 3L, diag)),
                covariances = covariance)
  n <- 100000L
  profile <- rep(1:2, each = n)
  set.seed(5)
  drawn <- .lta_draw_continuous(block, profile, c("a", "b", "c"))
  invisible(lapply(1:2, function(k) {
    rows <- profile == k
    # Sampling SD of a covariance entry is at most sqrt(9 / n) ~ 0.01 here.
    expect_lt(max(abs(stats::cov(drawn[rows, ]) - covariance[, , k])), 0.05)
    expect_lt(max(abs(colMeans(drawn[rows, ]) - means[k, ])), 0.03)
  }))
})

test_that("relabelling general logits preserves every implied probability", {
  set.seed(8)
  n_profiles <- 3L
  classes <- c("group_class_1", "group_class_2")
  terms <- c("(Intercept)", "z")
  profiles <- paste0("profile_", seq_len(n_profiles))
  pairs <- expand.grid(to = seq_len(n_profiles), from = seq_len(n_profiles))
  pairs <- pairs[pairs$to != pairs$from, ]
  triples <- expand.grid(to = seq_len(n_profiles), from = seq_len(n_profiles),
                         previous = seq_len(n_profiles))
  triples <- triples[triples$to != triples$from, ]
  names_of <- c(
    unlist(lapply(classes, function(h) outer(profiles[-n_profiles], terms, function(p, t)
      paste("initial", h, p, t, sep = ".")))),
    unlist(lapply(classes, function(h) outer(
      paste0(profiles[pairs$from], "->", profiles[pairs$to]), terms,
      function(p, t) paste("transition", h, p, t, sep = ".")))),
    unlist(lapply(classes, function(h) outer(
      paste0(profiles[triples$previous], "->", profiles[triples$from], "->",
             profiles[triples$to]), terms,
      function(p, t) paste("transition2", h, p, t, sep = ".")))))
  logits <- stats::setNames(stats::rnorm(length(names_of)), names_of)
  z <- 0.7
  value <- function(x, key) x[[paste0(key, ".(Intercept)")]] + z * x[[paste0(key, ".z")]]
  softmax <- function(v) exp(v) / sum(exp(v))
  initial <- function(x, h) softmax(c(vapply(profiles[-n_profiles], function(p)
    value(x, paste("initial", h, p, sep = ".")), numeric(1)), 0))
  move <- function(x, h, from) softmax(vapply(seq_len(n_profiles), function(to)
    if (to == from) 0 else value(x, paste("transition", h,
                                          paste0(profiles[from], "->", profiles[to]),
                                          sep = ".")), numeric(1)))
  move2 <- function(x, h, previous, from) softmax(vapply(seq_len(n_profiles), function(to)
    if (to == from) 0 else value(x, paste("transition2", h, paste0(
      profiles[previous], "->", profiles[from], "->", profiles[to]), sep = ".")),
    numeric(1)))
  order <- c(3L, 1L, 2L)
  class_order <- c(2L, 1L)
  moved <- .lta_relabel_classes(.lta_relabel_profiles(logits, order), class_order, classes)
  expect_setequal(names(moved), names(logits))
  invisible(lapply(seq_along(classes), function(i) {
    new_class <- classes[i]
    old_class <- classes[class_order[i]]
    expect_equal(unname(initial(moved, new_class)),
                 unname(initial(logits, old_class)[order]), tolerance = 1e-12)
    invisible(lapply(seq_len(n_profiles), function(j) {
      expect_equal(move(moved, new_class, j), move(logits, old_class, order[j])[order],
                   tolerance = 1e-12)
      lapply(seq_len(n_profiles), function(i_prev) {
        expect_equal(move2(moved, new_class, i_prev, j),
                     move2(logits, old_class, order[i_prev], order[j])[order],
                     tolerance = 1e-12)
      })
    }))
  }))
})

test_that("a permuted homogeneous fit aligns back to the original exactly", {
  skip_on_cran()
  data <- boot_lta_data()
  fit <- boot_lta_fit(data, n_group_classes = 2, model = "VVV")
  permuted <- .lta_permute_transitions(fit, c(2L, 1L), c(2L, 1L))
  expect_false(isTRUE(all.equal(.lta_bootstrap_estimates(permuted),
                                .lta_bootstrap_estimates(fit))))
  expect_equal(.lta_align_estimates(permuted, fit), .lta_bootstrap_estimates(fit),
               tolerance = 1e-12)
})

test_that("refits whose labels switched align to the original estimates", {
  skip_on_cran()
  data <- boot_lta_data()
  specs <- list(homogeneous = list(), covariance = list(transitions = "occasion",
                                                         model = "VVV"),
                second = list(transitions = "occasion", order = 2))
  invisible(lapply(specs, function(spec) {
    reference <- do.call(boot_lta_fit, c(list(data), spec))
    target <- .lta_bootstrap_estimates(reference)
    switched <- vapply(1:6, function(s) {
      set.seed(100 + s)
      persons <- sample(unique(data$id))
      shuffled <- data[order(match(data$id, persons), data$time), ]
      refit <- suppressWarnings(do.call(lta, c(list(shuffled, c("y1", "y2"), "id", 2,
                                                    time = "time", n_starts = 1, seed = s),
                                               spec)))
      expect_equal(refit$log_likelihood, reference$log_likelihood, tolerance = 1e-6)
      expect_equal(.lta_align_estimates(refit, reference)[names(target)], target,
                   tolerance = 1e-3)
      !identical(.multilpa_match_order(.lta_profile_signature(reference, reference),
                                       .lta_profile_signature(refit, reference)), 1:2)
    }, logical(1))
    # The check is only meaningful if some refit really came back relabelled.
    expect_true(any(switched))
  }))
})

test_that("transition bootstrap matches the Wald table's layout and scale", {
  skip_on_cran()
  data <- boot_lta_data()
  fits <- list(boot_lta_fit(data), boot_lta_fit(data, transitions = "occasion"))
  invisible(lapply(fits, function(fit) {
    wald <- parameter_inference(fit)
    boot <- parameter_inference(fit, data = data, method = "bootstrap", iter = 40,
                                n_starts = 1, seed = 3)
    expect_identical(names(boot), names(wald))
    expect_identical(nrow(boot), nrow(wald))
    expect_equal(boot$estimate, wald$estimate, tolerance = 1e-12)
    expect_true(all(is.na(boot$statistic) & is.na(boot$p_value) & is.na(boot$p_adjusted)))
    expect_identical(attr(boot, "method"), "bootstrap")
    expect_identical(attr(boot, "n_valid"), 40L)
    tested <- is.finite(wald$standard_error) & wald$standard_error > 0
    ratio <- boot$standard_error[tested] / wald$standard_error[tested]
    expect_true(all(ratio > 0.6 & ratio < 1.7))
    expect_true(all(boot$conf_low <= boot$conf_high))
  }))
})

test_that("the transition bootstrap is reproducible and leaves the caller's seed alone", {
  data <- boot_lta_data()
  fit <- boot_lta_fit(data)
  set.seed(99)
  before <- .Random.seed
  first <- parameter_inference(fit, data = data, method = "bootstrap", iter = 5,
                               n_starts = 1, seed = 4)
  expect_identical(.Random.seed, before)
  second <- parameter_inference(fit, data = data, method = "bootstrap", iter = 5,
                                n_starts = 1, seed = 4)
  expect_identical(first, second)
})

test_that("the bootstrap gives inference where Wald is refused", {
  skip_on_cran()
  data <- boot_lta_data()
  fit <- boot_lta_fit(data, transitions = "occasion", model = "VVV")
  expect_error(parameter_inference(fit), class = "latents_unsupported_inference")
  boot <- parameter_inference(fit, data = data, method = "bootstrap", iter = 30,
                              n_starts = 1, seed = 2)
  expect_true(all(is.finite(boot$standard_error) & boot$standard_error > 0))
  expect_identical(unique(boot$block), c("transition", "initial"))
})

test_that("weighted transition fits bootstrap with their weights", {
  skip_on_cran()
  data <- boot_lta_data()
  set.seed(6)
  data$w <- rep(sample(1:3, 120L, replace = TRUE), each = 4L)
  fit <- boot_lta_fit(data, weights = "w")
  boot <- parameter_inference(fit, data = data, method = "bootstrap", iter = 20,
                              n_starts = 1, seed = 1)
  expect_true(all(is.finite(boot$standard_error)))
  expect_equal(boot$estimate, parameter_inference(fit)$estimate, tolerance = 1e-12)
})

test_that("the transition bootstrap refuses bad requests by class", {
  data <- boot_lta_data()
  fit <- boot_lta_fit(data)
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 5),
               class = "latents_bad_argument")
  expect_error(parameter_inference(fit, data = data[-1L, ], method = "bootstrap",
                                   iter = 5),
               class = "latents_bad_inference_data")
  expect_error(parameter_inference(fit, data = data, method = "bootstrap", iter = 1))
})
