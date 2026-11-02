skip_on_cran()

# Audit 2026-11-01, slice additive-compare: group-class families, the tree
# table and compare_models(). See validation/audit-2026-11-01/additive-compare.md.

audit_family_data <- function() {
  truth <- list(means = rbind(c(0, 5), c(4, 1)),
                between = rbind(c(1, 0.5), c(0.6, 1.2)),
                within = c(2, 1), weights = c(0.6, 0.4))
  sizes <- rep(c(3L, 5L, 4L, 6L), 15L)
  additive_draw(truth, sizes, seed = 3)
}

audit_family_grid <- function() {
  suppressWarnings(enumerate_classes(audit_family_data(), c("y1", "y2"), "group",
                                     family = c("additive", "dispersion"),
                                     n_group_classes = 1:3, n_starts = 2, seed = 1))
}

test_that("a family grid refuses an unknown table or criterion with a class", {
  grid <- audit_family_grid()
  expect_error(get_results(grid, "bogus"), class = "latents_bad_argument")
  expect_error(get_results(grid, "candidates", criterion = "xx"),
               class = "latents_bad_argument")
  # Valid and partially matched names still work as before.
  expect_identical(get_results(grid, "cand"), get_results(grid, "candidates"))
  expect_identical(nrow(get_results(grid, "best", criterion = "aic")), 1L)
})

test_that("a family grid refuses invalid group-class counts before fitting", {
  data <- audit_family_data()
  bad <- list(c(0, 2), 2.5, NA, "2", integer())
  invisible(lapply(bad, function(count) {
    expect_error(enumerate_classes(data, c("y1", "y2"), "group", family = "additive",
                                   n_group_classes = count, n_starts = 1, seed = 1),
                 class = "latents_bad_argument")
  }))
})

test_that("compare_models() refuses a bad `cutoff` with a class", {
  shared <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
                model = "EEI", n_starts = 2, seed = 1)
  varying <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
                 model = "VVI", n_starts = 2, seed = 1)
  expect_error(compare_models(a = shared, b = varying, cutoff = 2),
               class = "latents_bad_argument")
  expect_error(compare_models(a = shared, b = varying, cutoff = c(0.5, 0.6)),
               class = "latents_bad_argument")
  # The LRT by hand: twice the log-likelihood difference on the parameter
  # difference (VVI adds two variances to EEI with two indicators).
  tested <- as.data.frame(compare_models(a = shared, b = varying, test = "lrt"))
  statistic <- 2 * (varying$log_likelihood - shared$log_likelihood)
  expect_equal(tested$lr_statistic, c(NA, statistic), tolerance = 1e-12)
  expect_identical(tested$lr_df, c(NA, 2L))
  expect_equal(tested$lr_p_value,
               c(NA, stats::pchisq(statistic, 2, lower.tail = FALSE)), tolerance = 1e-12)
})

test_that("group-class family fits get their own refusal from compare_models()", {
  data <- audit_family_data()
  first <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                    family = "additive", n_starts = 2, seed = 1)
  second <- suppressWarnings(multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                                      family = "additive", between_variance = "equal",
                                      n_starts = 2, seed = 1))
  # The class was always latents_bad_argument; the message used to say
  # compare_models() takes fits from multilpa(), which these are.
  expect_error(compare_models(a = second, b = first), class = "latents_bad_argument",
               regexp = "bootstrap_lrt")
  expect_error(compare_models(a = second, b = first, test = "lrt"),
               class = "latents_bad_argument", regexp = "bootstrap_lrt")
})

test_that("a tree's `from` and `to` take one count on a family grid", {
  grid <- audit_family_grid()
  expect_error(get_results(grid, "tree", from = c(1, 3), to = c(2, 3), model = "additive"),
               class = "latents_bad_argument")
  pair <- get_results(grid, "tree", from = 1, to = 3, model = "additive")
  expect_identical(unique(pair$from_candidate), "additive_groups1")
  expect_identical(unique(pair$to_candidate), "additive_groups3")
})

test_that("tree flows are a joint distribution with the shares as margins", {
  check_flows <- function(table) {
    pair <- paste(table$tree, table$from_candidate, table$to_candidate)
    expect_equal(as.vector(tapply(table$mass, pair, sum)),
                 rep(1, length(unique(pair))), tolerance = 1e-12)
    expect_equal(as.vector(tapply(table$share_of_from, paste(pair, table$from_class), sum)),
                 rep(1, length(unique(paste(pair, table$from_class)))), tolerance = 1e-12)
    expect_equal(as.vector(tapply(table$share_of_to, paste(pair, table$to_class), sum)),
                 rep(1, length(unique(paste(pair, table$to_class)))), tolerance = 1e-12)
    # The marginal of the flows is each row's share.
    from_mass <- tapply(table$mass, paste(pair, table$from_class), sum)
    from_share <- tapply(table$from_share, paste(pair, table$from_class), max)
    expect_equal(as.vector(from_mass), as.vector(from_share), tolerance = 1e-12)
    expect_true(all(table$mass >= 0))
  }
  family_grid <- audit_family_grid()
  check_flows(get_results(family_grid, "tree"))
  profile_grid <- enumerate_classes(subset(course_engagement, student <= 30),
                                    c("browse", "lectures"), "student",
                                    n_profiles = 1:3, n_group_classes = 1:2,
                                    model = "VVI", n_starts = 2, seed = 1)
  check_flows(get_results(profile_grid, "tree"))
  check_flows(get_results(profile_grid, "tree", level = "groups"))
  check_flows(get_results(profile_grid, "all")$tree)
})

test_that("additive fits reproduce the dense likelihood and their M-step means", {
  data <- audit_family_data()
  fit <- multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                  family = "additive", n_starts = 2, seed = 1)
  # Each group's ratings on indicator r: N(mu_hr 1, W_hr I + T_hr 11'), mixed
  # over the group classes; a dense Cholesky evaluation, independent of the
  # sufficient-statistic formula.
  groups <- unique(data$group)
  dense <- sum(vapply(groups, function(group) {
    ratings <- as.matrix(data[data$group == group, c("y1", "y2")])
    n <- nrow(ratings)
    joint <- vapply(1:2, function(h) {
      log(fit$group_probabilities[[h]]) + sum(vapply(1:2, function(r) {
        sigma <- diag(fit$within_variances[h, r], n) + fit$between_variances[h, r]
        root <- chol(sigma)
        z <- backsolve(root, ratings[, r] - fit$means[h, r], transpose = TRUE)
        -0.5 * n * log(2 * pi) - sum(log(diag(root))) - 0.5 * sum(z^2)
      }, numeric(1)))
    }, numeric(1))
    max(joint) + log(sum(exp(joint - max(joint))))
  }, numeric(1)))
  expect_equal(fit$log_likelihood, dense, tolerance = 1e-12)
  # Class mean = posterior average of the groups' conditional intercepts
  # (the M-step's fixed point), and a group's intercept mean mixes its
  # conditional means over its class posterior.
  posterior <- fit$group_posteriors
  conditional <- fit$conditional_intercept_means
  fixed_point <- t(vapply(1:2, function(h) {
    colSums(posterior[, h] * conditional[, , h]) / sum(posterior[, h])
  }, numeric(2)))
  expect_equal(unname(fixed_point), unname(fit$means), tolerance = 1e-8)
  expect_equal(unname(fit$intercept_means),
               unname(posterior[, 1] * conditional[, , 1] + posterior[, 2] * conditional[, , 2]),
               tolerance = 1e-12)
  # Hand counts: means, between and within variances (H rows when they vary,
  # one when shared) times d, plus H - 1 weights; H = 3, d = 2.
  counts <- c(additive = 16L, additive_equal = 12L, dispersion = 12L,
              additive_dispersion = 20L, additive_dispersion_equal = 16L)
  structures <- list(.additive_structure("additive", "varying"),
                     .additive_structure("additive", "equal"),
                     .additive_structure("dispersion", "equal"),
                     .additive_structure("additive_dispersion", "varying"),
                     .additive_structure("additive_dispersion", "equal"))
  expect_identical(vapply(structures, function(s) .additive_count_parameters(3L, 2L, s),
                          integer(1)), unname(counts))
})

test_that(".additive_nested() follows the block ranks of the families", {
  data <- audit_family_data()
  specs <- list(c("additive", "varying"), c("additive", "equal"),
                c("dispersion", "equal"), c("additive_dispersion", "varying"),
                c("additive_dispersion", "equal"))
  fits <- lapply(specs, function(spec) {
    suppressWarnings(multilpa(data, c("y1", "y2"), "group", n_group_classes = 2,
                              family = spec[1], between_variance = spec[2],
                              n_starts = 1, seed = 1))
  })
  nested <- outer(seq_along(fits), seq_along(fits), Vectorize(function(i, k) {
    .additive_nested(fits[[i]], fits[[k]])
  }))
  # Rows null, columns alternative, in the order of `specs`: a block shared
  # by every class is a special case of a class-specific one.
  expected <- rbind(c(TRUE, FALSE, FALSE, TRUE, FALSE),
                    c(TRUE, TRUE, FALSE, TRUE, TRUE),
                    c(FALSE, FALSE, TRUE, TRUE, TRUE),
                    c(FALSE, FALSE, FALSE, TRUE, FALSE),
                    c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_identical(nested, expected)
})

test_that("compare_models() reports the fits' own information criteria", {
  grid <- enumerate_classes(subset(course_engagement, student <= 30),
                            c("browse", "lectures"), "student",
                            n_profiles = 2:3, n_group_classes = 1:2,
                            model = "VVI", n_starts = 2, seed = 1)
  table <- as.data.frame(compare_models(grid))
  criteria <- do.call(rbind, lapply(table$model, function(name) {
    get_results(get_candidate(grid, name), "information_criteria", format = "wide")
  }))
  shared <- intersect(names(criteria), names(table))
  expect_gt(length(shared), 10L)
  expect_equal(table[shared], criteria[shared], tolerance = 1e-12, ignore_attr = TRUE)
  candidates <- get_results(grid, "candidates")
  matched <- candidates[match(table$model, candidates$candidate), shared]
  expect_equal(table[shared], matched, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(table$delta_bic, table$bic_groups - min(table$bic_groups))
  expect_equal(sum(table$bic_weight), 1)
})

test_that("additive inference refuses a bad `step` with a class", {
  fit <- suppressWarnings(multilpa(audit_family_data(), c("y1", "y2"), "group",
                                   n_group_classes = 2, family = "dispersion",
                                   n_starts = 1, seed = 1))
  invisible(lapply(list(-1, 0, NA_real_, Inf, c(1e-4, 1e-3), "a"), function(step) {
    expect_error(parameter_inference(fit, step = step), class = "latents_bad_argument")
  }))
})
