# Review 2026-10-26, INF-4 and INF-5: a label holding `: , [ ] %` -- an
# ordinal level, a continuous indicator name, a category -- must keep its
# standard errors, and every table must show the label's own text. A term that
# pairs two labels writes `:` and `%` inside each as `%3A` and `%25`, so the
# pair still splits at its one bare `:`.

# More than 1.5 s under CRAN's settings: skipped there to keep the check within
# its time limit; it runs in CI and locally.
skip_on_cran()

label_fit_data <- function() {
  withr::local_seed(1)
  data <- data.frame(g = rep(seq_len(20), each = 5),
                     a = c(stats::rnorm(50, -1.5), stats::rnorm(50, 1.5)),
                     b = stats::rnorm(100))
  names(data)[2:3] <- c("x,1", "y:2")
  data
}

ordinal_fit_data <- function() {
  withr::local_seed(12)
  n <- 160
  profile <- sample.int(2, n, TRUE)
  data <- data.frame(g = rep(seq_len(40), each = 4),
                     y = c(-1.5, 1.5)[profile] + stats::rnorm(n),
                     w = stats::rnorm(n))
  code <- pmin(3L, pmax(1L, as.integer(round(c(1.5, 2.5)[profile] +
                                                 stats::rnorm(n, 0, 0.6)))))
  data$o <- factor(c("low", "mid, ok", "high:[x]")[code],
                   levels = c("low", "mid, ok", "high:[x]"), ordered = TRUE)
  data
}

categorical_fit_data <- function() {
  withr::local_seed(5)
  n_groups <- 40
  data <- data.frame(g = rep(seq_len(n_groups), each = 5))
  group_class <- rep(stats::rbinom(n_groups, 1, 0.5), each = 5)
  profile <- stats::rbinom(nrow(data), 1, ifelse(group_class == 1, 0.8, 0.25))
  data$a <- stats::rnorm(nrow(data), ifelse(profile == 1, 2, 0), 1)
  data$b <- 0.5 * data$a + stats::rnorm(nrow(data), ifelse(profile == 1, -1, 0.5), 1)
  data$c <- factor(ifelse(stats::runif(nrow(data)) < ifelse(profile == 1, 0.7, 0.3),
                          "x, y", ifelse(stats::runif(nrow(data)) < 0.5, "p:q", "z")))
  data$w <- stats::rnorm(nrow(data))
  names(data)[names(data) == "a"] <- "a,1"
  names(data)[names(data) == "b"] <- "b:2"
  data
}

test_that("label escapes round-trip and a paired term splits at one colon", {
  labels <- c("plain", "a:b", "x, y", "[x]", "100%", "a%3Ab", "%25", "\r", "",
              "a,b:c[d]%e")
  expect_identical(.latents_parameter_decode(.latents_parameter_component(labels)),
                   labels)
  # Only the term separator and the escape character are escaped in a term.
  expect_identical(.latents_term_component(labels),
                   c("plain", "a%3Ab", "x, y", "[x]", "100%25", "a%253Ab", "%2525",
                     "\r", "", "a,b%3Ac[d]%25e"))
  terms <- .latents_term(labels, rev(labels))
  pieces <- strsplit(terms, ":", fixed = TRUE)
  expect_true(all(lengths(pieces[nzchar(rev(labels))]) == 2L))
  # Two pairs that join to the same text stay distinct.
  expect_false(identical(.latents_term("a:b", "c"), .latents_term("a", "b:c")))
})

test_that("continuous names with a comma or colon keep their errors and text", {
  data <- label_fit_data()
  fit <- quietly(multilpa(data, c("x,1", "y:2"), "g", n_profiles = 2,
                          n_group_classes = 1, n_starts = 2, seed = 1,
                          covariance_model = "full"))
  inference <- parameter_inference(fit)
  means <- subset(inference, parameter == "mean")
  expect_setequal(means$term, c("x,1", "y:2"))
  covariances <- subset(inference, parameter == "covariance")
  expect_setequal(covariances$term, c("x,1:x,1", "y%3A2:x,1", "y%3A2:y%3A2"))
  expect_true(all(is.finite(inference$standard_error)))
  # A covariance on the diagonal is a variance: an interval, no test.
  diagonal <- subset(covariances, term %in% c("x,1:x,1", "y%3A2:y%3A2"))
  expect_identical(nrow(diagonal), 4L)
  expect_true(all(is.na(diagonal$p_value)))
  expect_true(all(diagonal$conf_low > 0))
  profiles <- get_results(fit, "profiles")
  expect_setequal(profiles$indicator, c("x,1", "y:2"))
  expect_true(all(is.finite(profiles$mean_standard_error)))
  expect_true(all(is.finite(profiles$variance_standard_error)))
  expect_identical(anyDuplicated(names(coef(fit))), 0L)
  expect_identical(.multilpa_parameter_names(inference), names(coef(fit)))
  expect_identical(rownames(confint(fit)), names(coef(fit)))
})

test_that("ordinal levels with a comma or colon keep their text in every family", {
  data <- ordinal_fit_data()
  plain <- quietly(multilpa(data, c("y", "o"), "g", n_profiles = 2,
                            n_group_classes = 1, ordinal = "o", n_starts = 2,
                            seed = 1))
  covariate <- quietly(multilpa(data, c("y", "o"), "g", n_profiles = 2,
                                n_group_classes = 1, ordinal = "o",
                                profile_covariates = "w", n_starts = 2, seed = 1))
  expected <- c("o:mid, ok", "o:high%3A[x]")
  lapply(list(plain, covariate), function(fit) {
    inference <- parameter_inference(fit)
    intercepts <- subset(inference, parameter == "ordinal_intercept")
    expect_identical(intercepts$term, expected)
    expect_true(all(is.finite(intercepts$standard_error)))
    expect_true(all(paste0("measurement.ordinal_intercept.shared.", expected) %in%
                      names(coef(fit))))
    expect_identical(anyDuplicated(names(coef(fit))), 0L)
    ordinal <- get_results(fit, "ordinal")
    expect_setequal(ordinal$category, c("low", "mid, ok", "high:[x]"))
    located <- subset(ordinal, location != 0)
    expect_true(all(is.finite(located$location_standard_error)))
  })
})

test_that("categories and names are spelled one way by every family and accessor", {
  data <- categorical_fit_data()
  vars <- c("a,1", "b:2", "c")
  plain <- quietly(multilpa(data, vars, "g", n_profiles = 2, n_group_classes = 1,
                            categorical = "c", covariance_model = "full",
                            n_starts = 3, seed = 1, tol = 1e-10))
  covariate <- quietly(multilpa(data, vars, "g", n_profiles = 2, n_group_classes = 1,
                                categorical = "c", covariance_model = "full",
                                profile_covariates = "w", n_starts = 3, seed = 1,
                                tol = 1e-10))
  measurement_terms <- function(fit) {
    inference <- parameter_inference(fit)
    subset(inference, level == "measurement")
  }
  plain_terms <- measurement_terms(plain)
  covariate_terms <- measurement_terms(covariate)
  expect_identical(plain_terms$term, covariate_terms$term)
  expect_identical(plain_terms$parameter, covariate_terms$parameter)
  responses <- subset(plain_terms, parameter == "response")
  # The comma is the label's own text; the colon inside `p:q` is escaped
  # because the term pairs it with its indicator.
  expect_setequal(responses$term, c("c:x, y", "c:p%3Aq", "c:z"))
  expect_false(any(grepl("%2C", c(plain_terms$term, names(coef(plain)),
                                  names(coef(covariate))), fixed = TRUE)))
  lapply(list(plain, covariate), function(fit) {
    table <- get_results(fit, "responses")
    expect_setequal(table$category, c("x, y", "p:q", "z"))
    expect_true(all(is.finite(table$probability_standard_error)))
    profiles <- get_results(fit, "profiles")
    expect_true(all(is.finite(profiles$mean_standard_error)))
    expect_true(all(is.finite(profiles$variance_standard_error)))
    wanted <- "measurement.response.profile_1.c:x, y"
    expect_true(wanted %in% names(coef(fit)))
    interval <- confint(fit, parm = wanted)
    row <- subset(parameter_inference(fit),
                  parameter == "response" & outcome == "profile_1" & term == "c:x, y")
    expect_equal(unname(interval[1L, ]), c(row$conf_low, row$conf_high),
                 tolerance = 1e-12)
  })
  # The names a fit holds at a bound are spelled as coef() spells them.
  bounded <- plain
  bounded$min_probability <- 0.2
  held <- .multilpa_response_bound(bounded)
  expect_true(length(held$natural) > 0L)
  expect_true(all(held$natural %in% names(coef(plain))))
  expect_true(all(held$unconstrained %in% names(coef(plain, scale = "unconstrained"))))
  expect_true(any(grepl("c:x, y", held$natural, fixed = TRUE)))
})

# The same labels in the other families: a constrained covariance structure,
# a transition model, direct effects and a general latent transition model.

odd_continuous_data <- function() {
  withr::local_seed(1)
  data <- data.frame(g = rep(seq_len(30), each = 5),
                     a = c(stats::rnorm(75, -1.5), stats::rnorm(75, 1.5)),
                     w = stats::rnorm(150))
  data$b <- 0.4 * data$a + stats::rnorm(150) + 0.6 * data$w
  names(data)[match(c("a", "b", "w"), names(data))] <- c("x,1", "y:2", "w:z")
  data
}

odd_panel_data <- function() {
  withr::local_seed(2)
  waves <- 4L
  frame <- do.call(rbind, lapply(seq_len(80), function(person) {
    path <- Reduce(function(previous, wave) {
      if (stats::runif(1L) < 0.8) previous else 3L - previous
    }, seq_len(waves - 1L), init = sample(1:2, 1L), accumulate = TRUE)
    data.frame(person = person, wave = seq_len(waves), state = path)
  }))
  frame$a <- stats::rnorm(nrow(frame), c(0, 2.5)[frame$state])
  frame$b <- stats::rnorm(nrow(frame), c(0, 2)[frame$state]) + 0.3 * frame$a
  frame$q <- ifelse(stats::runif(nrow(frame)) < c(0.75, 0.25)[frame$state], "x, y",
                    ifelse(stats::runif(nrow(frame)) < 0.5, "p:q", "z"))
  frame$w <- stats::rnorm(nrow(frame))
  frame$state <- NULL
  names(frame)[match(c("a", "b"), names(frame))] <- c("x,1", "y:2")
  frame
}

test_that("a constrained structure's chart labels keep a comma or colon intact", {
  data <- odd_continuous_data()
  fit <- quietly(multilpa(data, c("x,1", "y:2"), "g", n_profiles = 2,
                          n_group_classes = 1, n_starts = 2, seed = 1,
                          volume = "varying", shape = "equal", orientation = "equal"))
  expect_identical(fit$covariance_structure, "VEE")
  chart <- names(coef(fit, scale = "unconstrained"))
  expect_true(all(c("measurement.log_shape_cholesky.shared.x,1:x,1",
                    "measurement.shape_cholesky.shared.y%3A2:x,1") %in% chart))
  expect_identical(anyDuplicated(chart), 0L)
  profiles <- get_results(fit, "profiles")
  expect_true(all(is.finite(profiles$mean_standard_error)))
  expect_true(all(is.finite(profiles$variance_standard_error)))
})

test_that("a transition model treats a diagonal covariance as a variance", {
  fit <- quietly(lta(odd_panel_data(), c("x,1", "y:2"), "person", 2L, time = "wave",
                     n_starts = 3, seed = 1, tol = 1e-9, covariance_model = "full"))
  expect_s3_class(fit, "multilpa_transitions")
  inference <- parameter_inference(fit)
  diagonal <- subset(inference, parameter == "covariance" &
                       term %in% c("x,1:x,1", "y%3A2:y%3A2"))
  expect_identical(nrow(diagonal), 4L)
  # No test of a variance against zero, and an interval that stays positive.
  expect_true(all(is.na(diagonal$p_value)))
  expect_true(all(diagonal$conf_low > 0))
  expect_setequal(subset(inference, parameter == "mean")$term, c("x,1", "y:2"))
  profiles <- get_results(fit, "profiles")
  expect_true(all(is.finite(profiles$variance_standard_error)))
})

test_that("direct effects of names with a comma or colon keep their errors", {
  data <- odd_continuous_data()
  fit <- quietly(multilpa(data, c("x,1", "y:2"), "g", n_profiles = 2,
                          n_group_classes = 1, n_starts = 2, seed = 1,
                          direct_effects = c("y:2 ~ w:z", "x,1 ~ w:z")))
  expect_true(all(c("measurement.direct_effect.shared.y%3A2:w%3Az",
                    "measurement.direct_effect.shared.x,1:w%3Az") %in%
                    names(coef(fit))))
  effects <- get_results(fit, "direct_effects")
  expect_identical(effects$indicator, c("y:2", "x,1"))
  expect_identical(effects$covariate, c("w:z", "w:z"))
  expect_true(all(is.finite(effects$standard_error)))
})

test_that("a latent transition model with covariates names labels by the same rule", {
  fit <- quietly(lta(odd_panel_data(), c("x,1", "y:2"), "person", 2L, time = "wave",
                     n_starts = 3, seed = 1, tol = 1e-9, covariance_model = "full",
                     transition_covariates = "w"))
  expect_s3_class(fit, "multilpa_lta")
  coefficients <- names(coef(fit))
  expect_true(all(c("mean.all.profile_1.x,1", "mean.all.profile_1.y:2",
                    "log_cholesky.all.profile_1.x,1:x,1",
                    "cholesky.all.profile_1.y%3A2:x,1",
                    "log_cholesky.all.profile_1.y%3A2:y%3A2") %in% coefficients))
  expect_identical(anyDuplicated(coefficients), 0L)
  expect_identical(rownames(vcov(fit)), coefficients)
})

test_that("a transition model spells its categories as the profile models do", {
  fit <- quietly(lta(odd_panel_data(), c("x,1", "q"), "person", 2L, time = "wave",
                     categorical = "q", n_starts = 3, seed = 1, tol = 1e-9))
  responses <- subset(parameter_inference(fit), parameter == "response")
  expect_setequal(responses$term, c("q:x, y", "q:p%3Aq", "q:z"))
  expect_true(all(c("measurement.response.profile_1.q:x, y",
                    "measurement.response.profile_1.q:p%3Aq") %in% names(coef(fit))))
  table <- get_results(fit, "responses")
  expect_true(all(is.finite(table$probability_standard_error)))
})
