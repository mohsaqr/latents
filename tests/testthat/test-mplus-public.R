test_that("single-level equal and varying LPA reproduce published Mplus examples", {
  invisible(lapply(c("7.9", "7.10"), function(example) {
    stopifnot(is.character(example))
    fixture <- readRDS(test_path("..", "fixtures", "mplus", paste0("public-", example, ".rds")))
    expect_equal(dim(fixture$data), c(500L, 6L))
    expect_false(anyNA(fixture$data))
    expect_identical(anyDuplicated(fixture$data), 0L)
    expect_equal(length(fixture$provenance$md5), 3L)
    expect_true(all(grepl("^https://www.statmodel.com/", fixture$provenance$urls)))
    fit <- multilpa(fixture$data, c("y1", "y2", "y3", "y4"), "cluster", 2L, 1L,
                     variance_model = fixture$variance_model, n_starts = 20L,
                     max_iter = 1000L, tol = 1e-12, seed = 912)
    expect_true(fit$converged)
    expect_false(fit$boundary)
    expect_equal(fit$n_failed_starts, 0L)
    expect_gte(fit$n_best_replicated, 2L)
    comparison <- compare_public_mplus(fit, fixture$expected)
    expect_true(all(comparison$passed), info = paste(capture.output(print(comparison)), collapse = "\n"))
    # A perturbed result must fail: the comparison is sensitive to regression.
    fit$means[, "y1"] <- fit$means[, "y1"] + 0.01
    expect_false(all(compare_public_mplus(fit, fixture$expected)$passed))
  }))
})
