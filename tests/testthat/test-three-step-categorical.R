# three_step() with a categorical distal outcome: class-specific category
# probabilities by the weighting methods and the ML step-three model, their
# intervals, pairwise differences and the overall Wald test. The comparison
# with Latent GOLD is in equivalence/test-latentgold-categorical-distal.R.

skip_on_cran()

categorical_case <- function(seed = 5, n = 900, centres = rbind(c(-2, 0), c(1.5, 2), c(1.5, -2))) {
  set.seed(seed)
  truth <- sample.int(3L, n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  shares <- rbind(c(0.7, 0.2, 0.1), c(0.2, 0.6, 0.2), c(0.1, 0.3, 0.6))
  data <- data.frame(y1 = centres[truth, 1L] + stats::rnorm(n),
                     y2 = centres[truth, 2L] + stats::rnorm(n))
  data$grade <- factor(vapply(truth, function(class) {
    sample(c("fail", "pass", "merit"), 1L, prob = shares[class, ])
  }, character(1)), levels = c("fail", "pass", "merit"))
  fit <- lpa(data, c("y1", "y2"), n_profiles = 3, model = "VVI", n_starts = 3, seed = 1)
  list(data = data, fit = fit)
}

test_that("BCH probabilities are the category shares under the BCH weights", {
  case <- categorical_case()
  result <- three_step(case$fit, case$data, "grade")
  expect_identical(nrow(result), 9L)
  expect_identical(result$category, rep(c("fail", "pass", "merit"), 3L))
  # The weights table, read back by unit and class, gives the same shares.
  weights <- get_results(case$fit, "bch_weights")
  weights$grade <- case$data$grade[weights$unit]
  shares <- aggregate(weight ~ class + grade, weights, sum)
  totals <- aggregate(weight ~ class, weights, sum)
  shares$probability <- shares$weight / totals$weight[match(shares$class, totals$class)]
  expected <- shares$probability[match(paste(result$class, result$category),
                                       paste(shares$class, shares$grade))]
  expect_equal(result$estimate, expected, tolerance = 1e-12)
  expect_equal(as.vector(tapply(result$estimate, result$class, sum)), rep(1, 3L),
               tolerance = 1e-12)
  expect_true(all(result$conf_low < result$estimate & result$estimate < result$conf_high))
  expect_true(all(result$conf_low > 0 & result$conf_high < 1))
})

test_that("the ML step-three model holds the class sizes at the average posteriors", {
  # Hard posteriors make the error matrix the identity: each unit belongs to
  # its assigned class, and the probabilities are the category shares.
  modal <- rep(1:2, c(30L, 20L))
  codes <- c(rep(1:3, c(15L, 10L, 5L)), rep(1:3, c(4L, 6L, 10L)))
  pieces <- list(posteriors = cbind(modal == 1L, modal == 2L) * 1, modal = modal,
                 n_classes = 2L)
  hard <- .multilpa_step_categorical_fit(pieces, outer(codes, 1:3, "==") * 1)
  expect_true(hard$converged)
  expect_equal(unname(hard$probabilities), rbind(c(15, 10, 5) / 30, c(4, 6, 10) / 20),
               tolerance = 1e-10)
  # On a real fit the sizes stay where step one put them, the probabilities
  # are close to BCH's and every row sums to one.
  case <- categorical_case(centres = rbind(c(-1.2, 0), c(1, 1.2), c(1, -1.2)))
  pieces <- .multilpa_level_assignments(case$fit, "individuals")
  soft <- .multilpa_step_categorical_fit(pieces,
                                         outer(as.integer(case$data$grade), 1:3, "==") * 1)
  expect_equal(soft$sizes, colMeans(pieces$posteriors), tolerance = 1e-12)
  bch <- three_step(case$fit, case$data, "grade")
  ml <- three_step(case$fit, case$data, "grade", method = "ml")
  expect_lt(max(abs(ml$estimate - bch$estimate) / bch$standard_error), 1)
  expect_equal(as.vector(tapply(ml$estimate, ml$class, sum)), rep(1, 3L), tolerance = 1e-12)
})

test_that("the ML information is the negative Hessian of its log likelihood", {
  skip_if_not_installed("numDeriv")
  case <- categorical_case(centres = rbind(c(-1.2, 0), c(1, 1.2), c(1, -1.2)))
  pieces <- .multilpa_level_assignments(case$fit, "individuals")
  indicators <- outer(as.integer(case$data$grade), seq_len(3L), "==") * 1
  fit <- .multilpa_step_categorical_fit(pieces, indicators)
  errors <- t(.multilpa_error_matrix(pieces)[, pieces$modal])
  softmax <- function(logits) exp(c(logits, 0)) / sum(exp(c(logits, 0)))
  log_likelihood <- function(parameters) {
    probabilities <- t(vapply(seq_len(3L), function(class) {
      softmax(parameters[(class - 1L) * 2L + 1:2])
    }, numeric(3)))
    sum(log(rowSums(errors * tcrossprod(indicators, probabilities) *
                      matrix(fit$sizes, nrow(indicators), 3L, byrow = TRUE))))
  }
  parameters <- as.vector(t(log(fit$probabilities[, 1:2] / fit$probabilities[, 3])))
  expect_lt(max(abs(numDeriv::grad(log_likelihood, parameters))), 1e-4)
  expect_equal(.multilpa_step_categorical_information(fit, indicators),
               -numDeriv::hessian(log_likelihood, parameters), tolerance = 1e-6)
})

test_that("results do not depend on how the categories are ordered or typed", {
  case <- categorical_case()
  reversed <- transform(case$data, grade = factor(grade, levels = c("merit", "pass", "fail")))
  original <- three_step(case$fit, case$data, "grade")
  turned <- three_step(case$fit, reversed, "grade")
  key <- function(table) paste(table$class, table$category)
  expect_equal(turned$estimate[match(key(original), key(turned))], original$estimate,
               tolerance = 1e-12)
  # The overall test does not depend on the baseline category.
  expect_equal(three_step(case$fit, reversed, "grade", contrast = "overall")$statistic,
               three_step(case$fit, case$data, "grade", contrast = "overall")$statistic,
               tolerance = 1e-8)
  as_text <- transform(case$data, grade = as.character(grade))
  expect_identical(sort(three_step(case$fit, as_text, "grade")$category, method = "radix"),
                   sort(original$category, method = "radix"))
  passed <- transform(case$data, passed = grade != "fail")
  logical_result <- three_step(case$fit, passed, "passed")
  expect_identical(unique(logical_result$category), c("FALSE", "TRUE"))
})

test_that("pairs and the overall test have one row per comparison", {
  case <- categorical_case()
  pairs <- three_step(case$fit, case$data, "grade", contrast = "pairs")
  expect_identical(nrow(pairs), 9L)
  expect_identical(attr(pairs, "adjust"), "BH")
  expect_equal(pairs$estimate[1L], with(three_step(case$fit, case$data, "grade"),
                                         estimate[4L] - estimate[1L]), tolerance = 1e-12)
  overall <- three_step(case$fit, case$data, "grade", contrast = "overall")
  expect_identical(overall$df, 4L)
  expect_lt(overall$p_value, 0.001)
  # A numeric outcome: with two classes the overall test is the squared z of
  # the one pair.
  two <- lpa(case$data, c("y1", "y2"), n_profiles = 2, model = "VVI", n_starts = 2, seed = 1)
  scored <- transform(case$data, score = y1 + stats::rnorm(nrow(case$data)))
  expect_equal(three_step(two, scored, "score", contrast = "overall")$statistic,
               three_step(two, scored, "score", contrast = "pairs")$statistic^2,
               tolerance = 1e-10)
  expect_identical(three_step(two, scored, "score", method = "ml", contrast = "overall")$df, 1L)
})

test_that("a BCH share outside [0, 1] is reported without an interval", {
  set.seed(10)
  truth <- sample.int(2L, 300L, replace = TRUE, prob = c(0.8, 0.2))
  data <- data.frame(y = stats::rnorm(300L, c(0, 1.2)[truth]),
                     z = stats::rnorm(300L, c(0, 1.2)[truth]))
  data$outcome <- factor(ifelse(stats::runif(300L) < c(0.02, 0.3)[truth], "rare", "common"))
  fit <- quietly(lpa(data, c("y", "z"), n_profiles = 2, n_starts = 2, seed = 1),
                 "latents_boundary")
  expect_warning(result <- three_step(fit, data, "outcome"), class = "latents_bch_out_of_range")
  outside <- subset(result, estimate < 0 | estimate > 1)
  expect_gt(nrow(outside), 0L)
  expect_true(all(is.na(outside$conf_low) & is.na(outside$conf_high)))
  expect_warning(overall <- three_step(fit, data, "outcome", contrast = "overall"),
                 class = "latents_bch_out_of_range")
  expect_true(is.na(overall$statistic))
  # The ML model keeps every probability in range; its maximum is then on
  # the bound, where it has no standard errors.
  expect_warning(ml <- three_step(fit, data, "outcome", method = "ml"),
                 class = "latents_boundary")
  expect_true(all(ml$estimate >= 0 & ml$estimate <= 1))
  expect_true(all(is.na(ml$standard_error)))
})

test_that("unusable categorical outcomes are refused", {
  case <- categorical_case()
  expect_error(three_step(case$fit, case$data, "grade", method = "ml", variance = "equal"),
               class = "latents_bad_argument")
  constant <- transform(case$data, grade = factor("pass"))
  expect_error(three_step(case$fit, constant, "grade"), class = "latents_bad_outcome")
  holes <- case$data
  holes$grade[3L] <- NA
  expect_error(three_step(case$fit, holes, "grade"), class = "latents_bad_argument")
})
