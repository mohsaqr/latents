skip_on_cran()

.fix1105_inference_case <- local({
  result <- NULL
  function() {
    if (!is.null(result)) return(result)
    withr::with_seed(15053, {
      n_groups <- 70L
      size <- 8L
      group <- rep(seq_len(n_groups), each = size)
      group_class <- rep(1:2, length.out = n_groups)
      profile <- 1L + stats::rbinom(length(group), 1L,
                                    c(0.3, 0.7)[group_class[group]])
      data <- data.frame(
        g = group, a = stats::rnorm(length(group), c(-0.8, 0.8)[profile]),
        b = stats::rnorm(length(group), c(-0.7, 0.7)[profile]),
        c = stats::rnorm(length(group), c(0.7, -0.7)[profile]),
        y = stats::rnorm(length(group), c(0, 1)[profile]),
        z = stats::rnorm(length(group)))
      staged <- quietly(fit_staged(data, c("a", "b", "c"), "g", 2L, 2L,
                                     n_starts = 3L, seed = 1L))
      result <<- list(data = data, staged = staged)
    })
    result
  }
})

test_that("a staged or held measurement cannot claim estimated step-one uncertainty", {
  case <- .fix1105_inference_case()
  conditional <- three_step(case$staged, case$data, "y", step_one = "fixed")
  expect_true(all(is.finite(conditional$standard_error)))
  invisible(lapply(c("bch", "ml", "modal"), function(method) {
    expect_error(three_step(case$staged, case$data, "y", method = method,
                            step_one = "estimated"),
                 class = "latents_unsupported_three_step")
  }))
  held <- quietly(multilpa(case$data, c("a", "b", "c"), "g", 2L, 1L,
                            start = starting_values(case$staged$stage_one),
                            fixed = "means", n_starts = 1L, seed = 1L))
  expect_false(isTRUE(held$staged))
  expect_error(three_step(held, case$data, "y", step_one = "estimated"),
               class = "latents_unsupported_three_step")
})

test_that("modal uncertainty does not require an unused step-one covariance", {
  case <- .fix1105_inference_case()
  fit <- case$staged$stage_one
  fixed <- three_step(fit, case$data, "y", method = "modal")
  estimated <- three_step(fit, case$data, "y", method = "modal", step_one = "estimated")
  expect_equal(estimated$estimate, fixed$estimate)
  expect_equal(estimated$standard_error, fixed$standard_error)
})

test_that("OPG probability and covariate SEs are invariant to grand transforms", {
  case <- .fix1105_inference_case()
  invisible(lapply(c(FALSE, TRUE), function(covariates) {
    arguments <- list(data = case$data, vars = c("a", "b", "c"), id = "g",
                       n_profiles = 2L, n_group_classes = 1L,
                       n_starts = 3L, seed = 1L, tol = 1e-11)
    if (covariates) arguments$profile_covariates <- "z"
    plain <- quietly(do.call(multilpa, arguments))
    reference <- parameter_inference(plain, vcov_type = "opg")
    keep <- reference$parameter %in% c("probability", "coefficient")
    invisible(lapply(list(list(centering = "grand"), list(scaling = "grand"),
                          list(centering = "grand", scaling = "grand")), function(transform) {
      fit <- quietly(do.call(multilpa, c(arguments, transform)))
      actual <- parameter_inference(fit, vcov_type = "opg")
      expect_equal(actual$estimate[keep], reference$estimate[keep], tolerance = 1e-5)
      expect_equal(actual$standard_error[keep], reference$standard_error[keep],
                    tolerance = 2e-4)
    }))
  }))
})

test_that("staged bootstrap checks the measurement sample before resampling", {
  case <- .fix1105_inference_case()
  expect_true(.multilpa_staged_same_data(case$staged))
  original <- case$data
  original$a <- original$a + 3
  measurement <- quietly(lpa(original, c("a", "b", "c"), 2L,
                               n_starts = 3L, seed = 1L))
  fit <- quietly(fit_staged(case$data, c("a", "b", "c"), "g", 2L, 2L,
                             measurement = measurement, n_starts = 3L, seed = 1L))
  expect_false(.multilpa_staged_same_data(fit))
  expect_error(parameter_inference(fit, method = "bootstrap", iter = 2L, seed = 1L),
               class = "latents_unsupported_inference")
  unknown <- case$staged
  unknown$stage_one$indicator_data <- NULL
  expect_false(.multilpa_staged_same_data(unknown))
  expect_error(parameter_inference(unknown, method = "bootstrap", iter = 2L),
               class = "latents_unsupported_inference")
})

test_that("same-data supplied measurement permits row order and identifier changes", {
  case <- .fix1105_inference_case()
  original <- case$data[rev(seq_len(nrow(case$data))), , drop = FALSE]
  original$renamed <- paste0("subject", original$g)
  measurement <- quietly(multilpa(original, c("a", "b", "c"), "renamed", 2L, 1L,
                                   n_starts = 3L, seed = 1L))
  fit <- quietly(fit_staged(case$data, c("a", "b", "c"), "g", 2L, 2L,
                             measurement = measurement, n_starts = 3L, seed = 1L))
  expect_true(.multilpa_staged_same_data(fit))
  boot <- quietly(parameter_inference(fit, method = "bootstrap", iter = 3L,
                                       n_starts = 2L, seed = 1L))
  expect_identical(attr(boot, "n_valid"), 3L)
  expect_identical(attr(boot, "stages_resampled"), 2L)
  expect_true(any(boot$parameter == "mean" & boot$standard_error > 0))
})

test_that("categorical sample provenance compares labels and preserves missing rows", {
  withr::with_seed(15056, {
    group <- rep(seq_len(40L), each = 6L)
    profile <- 1L + stats::rbinom(length(group), 1L, rep(c(.2, .8), each = 120L))
    data <- data.frame(g = group)
    vars <- paste0("v", 1:4)
    data[vars] <- lapply(vars, function(v) {
      factor(ifelse(stats::runif(length(group)) < c(.2, .8)[profile], "yes", "no"))
    })
    data$v1[c(2L, 8L)] <- NA
    other <- data[rev(seq_len(nrow(data))), , drop = FALSE]
    other[vars] <- lapply(other[vars], function(column) factor(column, levels = c("yes", "no")))
    measurement <- quietly(lca(other, vars, 2L, n_starts = 3L, seed = 1L, missing = "fiml"))
    fit <- quietly(fit_staged(data, vars, "g", 2L, 2L, categorical = vars,
                               measurement = measurement, n_starts = 3L, seed = 1L,
                               missing = "fiml"))
    expect_true(.multilpa_staged_same_data(fit))
    changed <- fit
    changed$stage_one$categorical_data[1L, "v2"] <-
      3L - changed$stage_one$categorical_data[1L, "v2"]
    expect_false(.multilpa_staged_same_data(changed))
    expect_error(parameter_inference(changed, method = "bootstrap", iter = 2L),
                 class = "latents_unsupported_inference")
  })
})

test_that("measurement provenance ignores grouping and sequence metadata", {
  case <- .fix1105_inference_case()
  data <- case$data
  data$time <- rep(seq_len(8L), 70L)
  original <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  original$person <- paste0("p", original$g)
  original$occasion <- original$time
  measurement <- quietly(multilpa(original, c("a", "b", "c"), "person", 2L, 1L,
                                   time = "occasion", n_starts = 1L, seed = 1L,
                                   max_iter = 20L))
  fit <- quietly(fit_staged(data, c("a", "b", "c"), "g", 2L, 2L,
                             time = "time", measurement = measurement,
                             n_starts = 1L, seed = 1L, max_iter = 20L))
  expect_true(.multilpa_staged_same_data(fit))
  changed <- fit
  changed$stage_one$time_values[1L] <- changed$stage_one$time_values[1L] + .5
  expect_true(.multilpa_staged_same_data(changed))
  changed <- fit
  indices <- which(changed$stage_one$time_values == 1L)[1:2]
  changed$stage_one$group_index[indices] <- rev(changed$stage_one$group_index[indices])
  expect_true(.multilpa_staged_same_data(changed))
})

test_that("step-three boundary checks use the whole design and unit-free coordinates", {
  withr::with_seed(15057, {
    z <- stats::rnorm(200L)
    design <- cbind(1, z)
    expect_error(.r3step_require_interior(matrix(c(0, 35), 2L), design),
                 class = "latents_boundary_fit")
    expect_error(.r3step_require_interior(matrix(c(-350000, .35), 2L),
                                          cbind(1, 1e6 + 100*z)),
                 class = "latents_boundary_fit")
    expect_error(.r3step_require_interior(matrix(c(-14.5, 0), 2L), design),
                 class = "latents_boundary_fit")
    design[1L, 2L] <- 20
    expect_no_error(.r3step_require_interior(matrix(c(0, 2), 2L), design))
    expect_no_error(.r3step_require_interior(matrix(c(-2e8, 2), 2L),
                                             cbind(1, 1e8 + z)))
  })
})

test_that("r3step refuses separated Wald tests under both covariance estimators", {
  withr::with_seed(15055, {
    n <- 400L
    profile <- 1L + stats::rbinom(n, 1L, .5)
    data <- data.frame(a = stats::rnorm(n, c(-2.5, 2.5)[profile]),
                        b = stats::rnorm(n, c(-2.5, 2.5)[profile]))
    fit <- quietly(lpa(data, c("a", "b"), 2L, n_starts = 3L, seed = 1L))
    data$x <- ifelse(fit$subject_profiles == 1L, 1, -1)*.6 + stats::rnorm(n, 0, .2)
    invisible(lapply(c("observed", "robust"), function(type) {
      expect_error(r3step(fit, data, "x", vcov_type = type), class = "latents_boundary_fit")
    }))
    data$x <- stats::rnorm(n)
    ordinary <- r3step(fit, data, "x")
    expect_true(all(is.finite(ordinary$standard_error)))
    data$x <- 1e8 + data$x
    shifted <- r3step(fit, data, "x")
    expect_equal(shifted$estimate[2L], ordinary$estimate[2L], tolerance = 1e-6)
    expect_equal(shifted$standard_error[2L], ordinary$standard_error[2L], tolerance = 1e-6)
  })
})

test_that("a rare transition regression cannot report finite sandwich Wald tests", {
  withr::with_seed(6, {
    n <- 80L
    n_occasions <- 4L
    n_profiles <- 3L
    states <- matrix(0L, n, n_occasions)
    states[, 1L] <- sample.int(n_profiles, n, TRUE)
    z <- matrix(stats::rnorm(n*n_occasions), n, n_occasions)
    invisible(lapply(2:n_occasions, function(t) {
      move <- stats::runif(n) < stats::plogis(-1.3 + .8*z[, t])
      states[, t] <<- ifelse(move, sample.int(n_profiles, n, TRUE), states[, t - 1L])
    }))
    means <- rbind(c(-1.2, -1.2), c(1.2, -1.2), c(0, 1.2))
    data <- data.frame(id = rep(seq_len(n), n_occasions),
                        time = rep(seq_len(n_occasions), each = n), z = as.vector(z),
                        a = stats::rnorm(n*n_occasions, means[as.vector(states), 1L]),
                        b = stats::rnorm(n*n_occasions, means[as.vector(states), 2L]))
    fit <- quietly(lta(data, c("a", "b"), id = "id", time = "time",
                        n_profiles = n_profiles, n_starts = 3L, seed = 1L))
    invisible(lapply(c("observed", "robust"), function(type) {
      expect_error(r3step(fit, data, transition_covariates = "z", vcov_type = type),
                   class = "latents_boundary_fit")
    }))
    occasion <- quietly(lta(data, c("a", "b"), id = "id", time = "time",
                             n_profiles = n_profiles, measurement = "occasion",
                             n_starts = 1L, seed = 1L, max_iter = 20L))
    expect_error(r3step(occasion, data, transition_covariates = "z"),
                 "invariant measurement", class = "latents_unsupported_three_step")
  })
})

test_that("binary category pairs count the same hypotheses as zero-one outcomes", {
  withr::with_seed(15058, {
    n <- 450L
    profile <- sample.int(3L, n, TRUE)
    data <- data.frame(a = stats::rnorm(n, c(-2, 0, 2)[profile]),
                        b = stats::rnorm(n, c(2, -2, 2)[profile]))
    data$binary <- factor(ifelse(stats::runif(n) < c(.35, .5, .65)[profile], "yes", "no"))
    data$numeric <- as.numeric(data$binary == "yes")
    fit <- quietly(lpa(data, c("a", "b"), 3L, n_starts = 3L, seed = 1L))
    invisible(lapply(c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY"), function(adjust) {
      categorical <- as.data.frame(three_step(fit, data, "binary", contrast = "pairs", adjust = adjust))
      numeric <- three_step(fit, data, "numeric", contrast = "pairs", adjust = adjust)
      chosen <- categorical$category == "yes"
      expect_equal(categorical$p_value[chosen], numeric$p_value, tolerance = 1e-10)
      expect_equal(categorical$p_adjusted[chosen], numeric$p_adjusted, tolerance = 1e-10)
    }))
    data$absent <- factor(ifelse(fit$subject_profiles == 1L, "yes", "no"))
    expect_warning(overall <- three_step(fit, data, "absent", method = "modal", contrast = "overall"),
                    class = "latents_boundary")
    expect_true(is.na(overall$p_value))
  })
})

test_that("weighted grand-transform inference uses the fitted weighted moments", {
  case <- .fix1105_inference_case()
  data <- case$data
  data$weight <- 1 + data$g %% 5L
  data$a[seq(3L, nrow(data), by = 17L)] <- NA_real_
  arguments <- list(data = data, vars = c("a", "b", "c"), id = "g",
                     n_profiles = 2L, n_group_classes = 1L,
                     n_starts = 3L, seed = 1L, tol = 1e-11,
                     weights = "weight", missing = "fiml")
  plain <- quietly(do.call(multilpa, arguments))
  grand <- quietly(do.call(multilpa, c(arguments,
                                       list(centering = "grand", scaling = "grand"))))
  p <- parameter_inference(plain)
  q <- parameter_inference(grand)
  keep <- p$parameter == "probability"
  expect_equal(q$standard_error[keep], p$standard_error[keep], tolerance = 1e-4)
  nuisance <- .multilpa_grand_transform(grand)
  expect_equal(colSums(nuisance$contributions), rep(0, length(nuisance$eta)),
                tolerance = 1e-8)
  raw <- as.matrix(data[c("a", "b", "c")])
  expected <- vapply(seq_len(ncol(raw)), function(j) {
    stats::weighted.mean(raw[, j], data$weight, na.rm = TRUE)
  }, numeric(1))
  expect_equal(nuisance$eta[seq_len(3L)], expected, tolerance = 1e-12)
})

test_that("zero Wald uncertainty cannot produce an infinite test statistic", {
  case <- .fix1105_inference_case()
  covariate <- quietly(multilpa(case$data, c("a", "b", "c"), "g", 2L, 1L,
                                 profile_covariates = "z", n_starts = 3L, seed = 1L))
  theta <- .multilpa_cov_encode(covariate)
  frame <- .multilpa_cov_inference_frame(covariate, theta,
                                          matrix(0, length(theta), length(theta)),
                                          .95, "observed", "none")
  expect_true(all(is.na(frame$statistic)))
  expect_true(all(is.na(frame$p_value)))
  one <- quietly(multilpa(case$data, c("a", "b", "c"), "g", 1L, 1L,
                           centering = "grand", n_starts = 1L, seed = 1L))
  result <- parameter_inference(one, vcov_type = "robust")
  zero <- result$standard_error == 0
  expect_true(any(zero))
  expect_true(all(is.na(result$statistic[zero])))
  expect_true(all(is.na(result$p_value[zero])))
})

test_that("covariate Wald inference warns about a material unfinished score", {
  withr::with_seed(5, {
    data <- do.call(rbind, lapply(seq_len(50L), function(g) {
      z <- stats::rnorm(6L)
      k <- 1L + stats::rbinom(6L, 1L, stats::plogis(.8*z))
      data.frame(g = g, z = z, y1 = stats::rnorm(6L, c(-1, 1)[k]),
                   y2 = stats::rnorm(6L, c(1, -1)[k]))
    }))
    loose <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
                                profile_covariates = "z", n_starts = 2L, seed = 1L,
                                tol = 1e-2))
    expect_warning(parameter_inference(loose), class = "latents_unconverged")
    tight <- quietly(multilpa(data, c("y1", "y2"), "g", 2L, 1L,
                                profile_covariates = "z", n_starts = 2L, seed = 1L,
                                tol = 1e-10))
    expect_no_warning(parameter_inference(tight))
  })
})

test_that("group-class r3step slopes survive a large covariate origin", {
  case <- .fix1105_inference_case()
  data <- case$data
  joint <- quietly(multilpa(data, c("a", "b", "c"), "g", 2L, 2L,
                             n_starts = 3L, seed = 1L))
  ordinary <- r3step(joint, data, "z", by_group_class = TRUE)
  data$z <- 1e8 + data$z
  shifted <- r3step(joint, data, "z", by_group_class = TRUE)
  keep <- ordinary$term == "z"
  expect_equal(shifted$estimate[keep], ordinary$estimate[keep], tolerance = 1e-5)
  expect_equal(shifted$standard_error[keep], ordinary$standard_error[keep], tolerance = 1e-5)
  expect_equal(attr(shifted, "log_likelihood"), attr(ordinary, "log_likelihood"),
                tolerance = 1e-8)
})
