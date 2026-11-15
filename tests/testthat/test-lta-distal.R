# three_step(occasion =): a person-level distal outcome of an lta() chain,
# related to the profile at one occasion or to the pair at two. The
# comparisons with independent references (a sum over every path, optim()
# and numDeriv) are in tests/equivalence/test-lta-distal.R; here, the
# structure of the result, identities and the refusals.

skip_on_cran()

chain_fixture <- function(n = 240L, occasions = 4L, seed = 7L) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(seed)
  state <- matrix(NA_integer_, n, occasions)
  state[, 1L] <- ifelse(stats::runif(n) < 0.5, 1L, 2L)
  for (t in seq_len(occasions)[-1L]) {  # simulation only
    moves <- stats::runif(n) < 0.2
    state[, t] <- ifelse(moves, 3L - state[, t - 1L], state[, t - 1L])
  }
  data <- data.frame(id = rep(seq_len(n), each = occasions),
                     time = rep(seq_len(occasions), n))
  s <- as.vector(t(state))
  data$y1 <- stats::rnorm(length(s), c(-1.2, 1.2)[s], 1)
  data$y2 <- stats::rnorm(length(s), c(-1, 1)[s], 1)
  z <- stats::rnorm(n, c(0, 2)[state[, 3L]], 1)
  data$z <- rep(z, each = occasions)
  data$group <- rep(ifelse(z > 1, "high", "low"), each = occasions)
  data$varying <- stats::rnorm(nrow(data))
  data
}

chain_fit <- function(data, ...) {
  lta(data, c("y1", "y2"), "id", n_profiles = 2, time = "time", n_starts = 2,
      seed = 1, ...)
}

test_that("a distal outcome at one occasion and at a pair has the documented layout", {
  data <- chain_fixture()
  fit <- chain_fit(data)
  one <- three_step(fit, data, "z", occasion = 3)
  expect_identical(attr(one, "occasion"), 3L)
  expect_identical(unique(one$level), "persons")
  expect_identical(one$class, 1:2)
  expect_match(attr(one, "title"), "profile at occasion 3", fixed = TRUE)
  pair <- three_step(fit, data, "z", occasion = c(2, 3), method = "ml")
  expect_identical(attr(pair, "occasion"), c(2L, 3L))
  expect_identical(pair$class, 1:4)
  expect_identical(pair$from, c(1L, 1L, 2L, 2L))
  expect_identical(pair$to, c(1L, 2L, 1L, 2L))
  expect_true(all(c("variance", "effective_n") %in% names(pair)))
  expect_s3_class(pair, "latents_table")
  expect_output(print(pair), "From")
  contrasts <- three_step(fit, data, "z", occasion = c(2, 3), contrast = "pairs")
  expect_true(all(c("from", "to", "reference_from", "reference_to") %in% names(contrasts)))
  expect_identical(nrow(contrasts), 6L)
  overall <- three_step(fit, data, "z", occasion = 3, contrast = "overall")
  expect_identical(overall$df, 1L)
  categorical <- three_step(fit, data, "group", occasion = 3)
  expect_identical(categorical$category, rep(c("high", "low"), 2L))
})

test_that("the pair's modal means average back to the later occasion's", {
  data <- chain_fixture()
  fit <- chain_fit(data)
  # Under modal assignment the effective n of a cell is its count, so the
  # count-weighted pair means over the earlier profile are the later
  # occasion's means.
  modal_pair <- three_step(fit, data, "z", occasion = c(2, 3), method = "modal")
  modal_one <- three_step(fit, data, "z", occasion = 3, method = "modal")
  marginal <- vapply(1:2, function(l) {
    rows <- modal_pair$to == l
    sum(modal_pair$estimate[rows] * modal_pair$effective_n[rows]) /
      sum(modal_pair$effective_n[rows])
  }, numeric(1))
  expect_equal(marginal, modal_one$estimate, tolerance = 1e-12)
})

test_that("the class means move with the outcome's origin and units", {
  data <- chain_fixture()
  fit <- chain_fit(data)
  shifted <- transform(data, z = 100 + 4 * z)
  for (method in c("bch", "ml")) {
    base <- three_step(fit, data, "z", occasion = 3, method = method)
    moved <- three_step(fit, shifted, "z", occasion = 3, method = method)
    expect_equal(moved$estimate, 100 + 4 * base$estimate, tolerance = 1e-6)
    expect_equal(moved$standard_error, 4 * base$standard_error, tolerance = 1e-5)
  }
})

test_that("the chain lets ML use a person whose chosen occasion is not measured", {
  data <- chain_fixture()
  gapped <- data[!(data$id %in% 1:20 & data$time == 3L), ]
  fit <- chain_fit(gapped, occasions = "grid")
  ml <- three_step(fit, gapped, "z", occasion = 3, method = "ml")
  bch <- three_step(fit, gapped, "z", occasion = 3)
  expect_equal(sum(ml$effective_n), 240, tolerance = 1e-9)
  expect_lt(sum(bch$effective_n), 221)
})

test_that("a distal outcome of a chain refuses what it does not define", {
  data <- chain_fixture(n = 160L)
  fit <- chain_fit(data)
  expect_error(three_step(fit, data, "z", occasion = 0), class = "latents_bad_argument")
  expect_error(three_step(fit, data, "z", occasion = 2.5), class = "latents_bad_argument")
  expect_error(three_step(fit, data, "z", occasion = c(3, 2)), class = "latents_bad_argument")
  expect_error(three_step(fit, data, "z", occasion = c(1, 2, 3)),
               class = "latents_bad_argument")
  expect_error(three_step(fit, data, "z", occasion = 5), class = "latents_bad_argument")
  expect_error(three_step(fit, data, "z", occasion = "3"), class = "latents_bad_argument")
  expect_error(three_step(fit, data, "z", occasion = 3, level = "groups"),
               class = "latents_bad_argument")
  expect_error(three_step(fit, data, "varying", occasion = 3), class = "latents_bad_outcome")
  expect_error(three_step(fit, data, "y1", occasion = 3), class = "multilpa_indicator_reused")
  expect_error(three_step(fit, data, "time", occasion = 3), class = "latents_bad_outcome")
  expect_error(three_step(fit, data, "z", occasion = 3, step_one = "estimated"),
               class = "latents_unsupported_three_step")
  expect_error(three_step(fit, subset(data, id > 1), "z", occasion = 3),
               class = "latents_bad_inference_data")
  profiles <- multilpa(data, c("y1", "y2"), "id", n_profiles = 2, n_group_classes = 1,
                       n_starts = 1, seed = 1)
  expect_error(three_step(profiles, data, "z", occasion = 3), class = "latents_bad_argument")
  covariate_data <- transform(data, w = sin(seq_len(nrow(data))))
  with_covariates <- chain_fit(covariate_data, transition_covariates = "w")
  expect_error(three_step(with_covariates, covariate_data, "z", occasion = 3),
               class = "latents_unsupported_three_step")
  # Two group classes in so small a sample can put a probability on its bound;
  # the refusal is the point here.
  two <- suppressWarnings(chain_fit(data, n_group_classes = 2))
  expect_error(three_step(two, data, "z", occasion = 3),
               class = "latents_unsupported_three_step")
})
