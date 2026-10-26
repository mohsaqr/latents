# Joint items: local dependence between two categorical indicators
# (`local_dependence = "a:b"`). Invariances, the dependence table recomputed
# from the estimates, prediction, the verbs that refit, the documented error
# classes, and two in-memory sabotages that these checks must catch. The
# comparisons with Latent GOLD, Mplus and poLCA are in
# equivalence/test-features-joint-items.R.

joint_items_data <- function(seed = 11L, n = 300L) {
  if (exists(".Random.seed", globalenv())) {
    saved <- get(".Random.seed", globalenv())
    on.exit(assign(".Random.seed", saved, globalenv()), add = TRUE)
  }
  set.seed(seed)
  class <- sample(1:2, n, TRUE, prob = c(0.55, 0.45))
  a <- vapply(class, \(j) sample(1:3, 1L, prob = list(c(0.6, 0.3, 0.1), c(0.1, 0.3, 0.6))[[j]]),
              integer(1))
  # b follows a within each class, more strongly in class 1.
  keep <- stats::runif(n) < c(0.7, 0.4)[class]
  b <- ifelse(keep, as.integer(a >= 2), stats::rbinom(n, 1, c(0.3, 0.7)[class]))
  data.frame(a = a, b = b,
             c = stats::rbinom(n, 1, c(0.8, 0.25)[class]),
             e = stats::rbinom(n, 1, c(0.75, 0.3)[class]),
             f = stats::rbinom(n, 1, c(0.7, 0.2)[class]))
}

joint_items_fit <- function(data, vars = c("a", "b", "c", "e", "f"), pairs = "a:b", ...) {
  quietly(lca(data, vars, 2, local_dependence = pairs, n_starts = 6, seed = 1, tol = 1e-12,
              max_iter = 5000, ...))
}

test_that("row order, pair order and column order leave the joint model unchanged", {
  data <- joint_items_data()
  fit <- joint_items_fit(data)
  rows <- rev(seq_len(nrow(data)))
  permuted <- joint_items_fit(data[rows, ])
  expect_equal(permuted$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  expect_identical(permuted$n_parameters, fit$n_parameters)
  # The posteriors of a row do not depend on where it sits (up to labels).
  back <- order(rows)
  same <- abs(permuted$subject_posteriors[back, 1] - fit$subject_posteriors[, 1])
  swapped <- abs(permuted$subject_posteriors[back, 2] - fit$subject_posteriors[, 1])
  # Both fits stop at a relative change of 1e-12, a few 1e-6 apart in a posterior.
  expect_lt(min(max(same), max(swapped)), 1e-5)
  backward <- joint_items_fit(data, pairs = "b:a")
  expect_equal(backward$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  expect_identical(backward$n_parameters, fit$n_parameters)
  reordered <- joint_items_fit(data[c("f", "e", "c", "b", "a")], vars = c("f", "e", "b", "c", "a"))
  expect_equal(reordered$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  # A swapped pair reports the same joint probabilities, items swapped.
  table <- get_results(fit, "local_dependence")
  mirrored <- get_results(backward, "local_dependence")
  expect_setequal(round(table$probability, 6), round(mirrored$probability, 6))
  expect_identical(unique(mirrored$item_1), "b")
})

test_that("multilpa() with categorical items fits the same joint model as lca()", {
  data <- joint_items_data()
  vars <- c("a", "b", "c", "e", "f")
  through_lca <- joint_items_fit(data)
  through_multilpa <- quietly(multilpa(data, vars, NULL, 2, 1, categorical = vars,
                                       local_dependence = "a:b", n_starts = 6, seed = 1,
                                       tol = 1e-12, max_iter = 5000))
  expect_equal(through_multilpa$log_likelihood, through_lca$log_likelihood, tolerance = 1e-9)
  expect_equal(through_multilpa$n_parameters, through_lca$n_parameters)
  expect_identical(through_multilpa$vars, c("a:b", "c", "e", "f"))
})

test_that("relabelled categories and recoded numbers leave the likelihood unchanged", {
  data <- joint_items_data()
  fit <- joint_items_fit(data)
  relabelled <- transform(data, a = c("high", "low", "mid")[a], b = c("yes", "no")[b + 1L])
  relabelled_fit <- joint_items_fit(relabelled)
  expect_equal(relabelled_fit$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
  expect_identical(relabelled_fit$n_parameters, fit$n_parameters)
  expect_setequal(relabelled_fit$categorical_levels[["a:b"]],
                  unique(paste(relabelled$a, relabelled$b, sep = ":")))
  # Factors with permuted level order are the same data.
  factors <- transform(data, a = factor(a, levels = c(3, 1, 2)), b = factor(b, levels = 1:0))
  expect_equal(joint_items_fit(factors)$log_likelihood, fit$log_likelihood, tolerance = 1e-9)
})

# Number codes that collide when written next to each other without a
# separator: (1, 10) and (11, 0) both read "110", (1, 1) and (11, "") "11".
joint_items_numeric_ok <- function() {
  data <- joint_items_data()
  data$a <- c(1, 11, 11)[data$a]
  data$b <- ifelse(data$f == 1, c(0, 10)[data$b + 1L], 1)
  combinations <- unique(data.frame(a = data$a, b = data$b))
  fit <- joint_items_fit(data, vars = c("a", "b", "c", "e"))
  table <- get_results(fit, "local_dependence")
  identical(length(fit$categorical_levels[["a:b"]]), nrow(combinations)) &&
    identical(nrow(table), 2L * nrow(combinations)) &&
    setequal(paste(table$category_1, table$category_2),
             paste(combinations$a, combinations$b))
}

test_that("number codes stay distinct combinations", {
  expect_true(joint_items_numeric_ok())
})

# The dependence table recomputed from the joint estimates alone.
joint_items_table_ok <- function(fit) {
  table <- get_results(fit, "local_dependence")
  joint <- fit$response_probabilities[["a:b"]]
  labels <- strsplit(colnames(joint), ":", fixed = TRUE)
  first <- vapply(labels, `[`, "", 1L)
  second <- vapply(labels, `[`, "", 2L)
  expected <- do.call(rbind, lapply(seq_len(nrow(joint)), function(j) {
    p <- joint[j, ]
    independent <- vapply(seq_along(p), \(i) sum(p[first == first[i]]) * sum(p[second == second[i]]),
                          numeric(1))
    data.frame(key = paste(j, first, second), probability = unname(p),
               independent = independent)
  }))
  at <- match(paste(table$profile, table$category_1, table$category_2), expected$key)
  !anyNA(at) && nrow(table) == nrow(expected) &&
    isTRUE(all.equal(table$probability, expected$probability[at], tolerance = 1e-12)) &&
    isTRUE(all.equal(table$independent, expected$independent[at], tolerance = 1e-12)) &&
    isTRUE(all.equal(table$ratio, table$probability / table$independent, tolerance = 1e-12))
}

test_that("the dependence table is the joint table beside the product of its margins", {
  fit <- joint_items_fit(joint_items_data())
  expect_true(joint_items_table_ok(fit))
  table <- get_results(fit, "local_dependence")
  totals <- aggregate(cbind(probability, independent) ~ profile, table, sum)
  expect_equal(totals$probability, c(1, 1), tolerance = 1e-12)
  # The joint table is the fit's joint indicator, no more and no less.
  expect_identical(nrow(table), 2L * length(fit$categorical_levels[["a:b"]]))
})

test_that("integer weights are duplicated rows", {
  data <- joint_items_data(n = 200L)
  data$w <- rep(1:3, length.out = nrow(data))
  weighted <- joint_items_fit(data, weights = "w")
  rows <- data[rep(seq_len(nrow(data)), data$w), ]
  duplicated_fit <- joint_items_fit(rows)
  # latents rescales weights to sum to the number of rows.
  expect_equal(weighted$log_likelihood, duplicated_fit$log_likelihood * nrow(data) / nrow(rows),
               tolerance = 1e-8)
  expect_identical(weighted$n_parameters, duplicated_fit$n_parameters)
  order <- if (abs(weighted$profile_probabilities[1, 1] -
                   duplicated_fit$profile_probabilities[1, 1]) < 0.05) 1:2 else 2:1
  expect_equal(unname(weighted$response_probabilities[["a:b"]]),
               unname(duplicated_fit$response_probabilities[["a:b"]][order, ]), tolerance = 1e-5)
})

test_that("prediction joins the pair in new data", {
  data <- joint_items_data()
  fit <- joint_items_fit(data)
  expect_equal(predict(fit, newdata = data, type = "posterior"),
               predict(fit, type = "posterior"), tolerance = 1e-12)
  expect_equal(as.matrix(predict(fit, newdata = data, type = "posterior")[
    c("probability_profile_1", "probability_profile_2")]),
    unname(fit$subject_posteriors), tolerance = 1e-10, ignore_attr = TRUE)
  # A combination the fit never saw has no estimated probability.
  unseen <- data[1L, ]
  unseen$a <- 4L
  expect_error(predict(fit, unseen, type = "posterior"), class = "latents_bad_data")
  # A row observing one item of the pair contributes that item, the pair's
  # other item summed out (review 2026-10-26 round 2, C2-2). The old
  # expectation here -- a latents_partial_pair warning, and the row
  # predicting exactly like one with the pair missing -- was the bias itself.
  gappy <- data
  gappy$b[1:4] <- NA
  fiml <- NULL
  expect_no_warning(fiml <- joint_items_fit(gappy, missing = "fiml"),
                    class = "latents_partial_pair")
  partial <- data[1:2, ]
  partial$a[1] <- NA
  prediction <- predict(fiml, partial, type = "posterior")
  joint <- fiml$response_probabilities[["a:b"]]
  second <- sub("^[^:]*:", "", fiml$categorical_levels[["a:b"]])
  item <- \(name) fiml$response_probabilities[[name]][, as.character(partial[[name]][1])]
  shares <- colMeans(fiml$subject_posteriors)
  kernel <- shares * rowSums(joint[, second == as.character(partial$b[1]), drop = FALSE]) *
    item("c") * item("e") * item("f")
  expect_equal(unname(unlist(prediction[1, c("probability_profile_1", "probability_profile_2")])),
               unname(kernel / sum(kernel)), tolerance = 1e-10)
  without_pair <- quietly(predict(fiml, transform(partial, a = NA, b = NA), type = "posterior"))
  expect_false(isTRUE(all.equal(prediction[1, ], without_pair[1, ], check.attributes = FALSE)))
})

test_that("refits, enumeration, sensitivity and comparison keep the pair", {
  skip_on_cran()
  data <- joint_items_data()
  fit <- joint_items_fit(data)
  frame <- .multilpa_resolve_data(fit, NULL)
  refit <- quietly(do.call(multilpa, c(list(data = frame, n_starts = 6, seed = 1, tol = 1e-12),
                                       .multilpa_refit_arguments(fit))))
  expect_equal(refit$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  expect_equal(refit$n_parameters, fit$n_parameters)
  labels <- .latents_categorical_labels(fit)
  expect_identical(colnames(labels), c("a", "b", "c", "e", "f"))
  expect_identical(unname(labels[, "a"]), as.character(data$a))
  grid <- quietly(enumerate_lca(data, c("a", "b", "c", "e", "f"), 1:2, local_dependence = "a:b",
                                seed = 1, n_starts = 6, tol = 1e-12))
  candidate <- candidate_fit(grid, 2)
  expect_false(is.null(candidate$local_dependence))
  expect_equal(candidate$log_likelihood, fit$log_likelihood, tolerance = 1e-8)
  stable <- quietly(sensitivity(fit, seeds = 1:2))
  expect_equal(as.data.frame(stable)$log_likelihood, rep(fit$log_likelihood, 2),
               tolerance = 1e-8)
  plain <- quietly(lca(data, c("a", "b", "c", "e", "f"), 2, n_starts = 6, seed = 1))
  compared <- as.data.frame(compare_models(plain, fit))
  expect_identical(nrow(compared), 2L)
  expect_equal(compared$n_parameters, c(plain$n_parameters, fit$n_parameters))
  expect_equal(compared$log_likelihood, c(plain$log_likelihood, fit$log_likelihood))
  boot <- quietly(parameter_inference(fit, method = "bootstrap", iter = 5, seed = 1,
                                      n_starts = 3))
  terms <- as.data.frame(boot)$term
  # Each class's joint probabilities are bootstrapped as the pair's terms.
  # The joined indicator's colon is escaped within its component; the final
  # colon still separates that indicator from its response category.
  joint <- !is.na(terms) & startsWith(terms, "a%3Ab:")
  expect_identical(sum(joint),
                   2L * length(fit$categorical_levels[["a:b"]]))
  expect_identical(sum(joint), 12L)
  expect_identical(anyDuplicated(.multilpa_parameter_names(as.data.frame(boot))), 0L)
})

test_that("bootstrap inference survives a resample without a rare combination", {
  skip_on_cran()
  # Finding 2 (equivalence/features/joint-items/findings/): a resample that
  # misses a combination observed twice used to stop with an internal,
  # unclassed error. Such a replicate is now dropped and counted.
  data <- joint_items_data()
  data$b[data$a == 1] <- 0
  data$b[which(data$a == 1)[1:2]] <- 1
  fit <- joint_items_fit(data)
  dropped <- FALSE
  boot <- withCallingHandlers(
    quietly(parameter_inference(fit, method = "bootstrap", iter = 20, seed = 1, n_starts = 3)),
    latents_bootstrap_dropped = function(condition) {
      dropped <<- TRUE
      invokeRestart("muffleWarning")
    })
  expect_true(dropped)
  expect_true(all(is.finite(boot$standard_error)))
})

test_that("pairs that cannot be joined raise their documented classes", {
  data <- joint_items_data(n = 120L)
  fit_with <- function(pairs, frame = data, ...) {
    lca(frame, c("a", "b", "c", "e"), 2, local_dependence = pairs, n_starts = 1, seed = 1, ...)
  }
  expect_error(fit_with("a:a"), class = "latents_bad_argument")
  expect_error(fit_with(c("a:b", "b:c")), class = "latents_bad_argument")
  expect_error(fit_with("a:b:c"), class = "latents_bad_argument")
  expect_error(fit_with("a:"), class = "latents_bad_argument")
  expect_error(fit_with(NA_character_), class = "latents_bad_argument")
  expect_error(fit_with("a:z"), class = "latents_bad_argument")
  expect_error(lca(data, c("a", "a", "b", "c"), 2, local_dependence = "a:b", n_starts = 1),
               class = "latents_bad_data")
  colon <- transform(data, a = c("x:1", "y", "z")[a])
  expect_error(fit_with("a:b", colon), class = "latents_bad_data")
  collision <- data
  collision[["a:b"]] <- 1
  expect_error(fit_with("a:b", collision), class = "latents_bad_data")
  expect_error(fit_with("a:b", transform(data, b = 1L)), class = "latents_bad_categorical")
  continuous <- transform(data, score = stats::qnorm(seq_len(nrow(data)) / (nrow(data) + 1)))
  expect_error(multilpa(continuous, c("a", "score", "c"), NULL, 2, 1, categorical = c("a", "c"),
                        local_dependence = "a:score", n_starts = 1, seed = 1),
               class = "latents_bad_argument")
  expect_error(quietly(multilpa(data, c("a", "b", "c"), NULL, 2, 1, categorical = "c",
                                ordinal = c("a", "b"), local_dependence = "a:b",
                                n_starts = 1, seed = 1)),
               class = "latents_bad_argument")
  expect_error(quietly(multilpa(data, c("a", "b", "c"), NULL, 2, 1, categorical = "c",
                                count = c("a", "b"), local_dependence = "a:b",
                                n_starts = 1, seed = 1)),
               class = "latents_bad_argument")
  expect_error(get_results(quietly(lca(data, c("a", "b", "c"), 2, n_starts = 1, seed = 1)),
                           "local_dependence"),
               class = "latents_no_local_dependence")
})

test_that("in-memory sabotage of the join and of the table is caught", {
  skip_on_cran()
  replace_in_namespace <- function(name, value) {
    original <- get(name, asNamespace("latents"))
    utils::assignInNamespace(name, value, "latents")
    original
  }
  # A join that writes the two codes next to each other without a separator.
  join <- get(".latents_join_pairs", asNamespace("latents"))
  sabotaged_join <- function(data, pairs, ...) {
    out <- join(data, pairs, ...)
    invisible(lapply(seq_len(nrow(pairs)), function(i) {
      name <- pairs$name[i]
      if (all(c(pairs$item_1[i], pairs$item_2[i]) %in% names(data))) {
        glued <- ifelse(is.na(out[[name]]), NA_character_,
                        paste0(data[[pairs$item_1[i]]], data[[pairs$item_2[i]]]))
        attributes_kept <- attributes(out[[name]])
        out[[name]] <<- factor(glued)
        attr(out[[name]], "latents_pair") <<- attributes_kept$latents_pair
      }
    }))
    out
  }
  original <- replace_in_namespace(".latents_join_pairs", sabotaged_join)
  # A sabotaged join may also fail outright; either way the check does not pass.
  detected <- tryCatch(!joint_items_numeric_ok(), error = function(e) TRUE)
  utils::assignInNamespace(".latents_join_pairs", original, "latents")
  expect_true(detected)
  expect_true(joint_items_numeric_ok())
  # A table whose independence column uses the pooled margins of all classes.
  fit <- joint_items_fit(joint_items_data())
  table_function <- get(".latents_local_dependence_table", asNamespace("latents"))
  sabotaged_table <- function(x) {
    out <- table_function(x)
    out$independent <- stats::ave(out$independent, out$category_1, out$category_2)
    out$ratio <- out$probability / out$independent
    out
  }
  original <- replace_in_namespace(".latents_local_dependence_table", sabotaged_table)
  detected <- !joint_items_table_ok(fit)
  utils::assignInNamespace(".latents_local_dependence_table", original, "latents")
  expect_true(detected)
  expect_true(joint_items_table_ok(fit))
})
