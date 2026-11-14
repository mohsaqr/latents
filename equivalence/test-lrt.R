# compare_models(test = "lrt") against mclust (covariance structures), glca's
# gofglca() (a membership covariate) and Mplus (LTA measurement invariance).
# The references are stored in lrt/reference.rds by lrt/make-kit.R. One block
# per reference program, by the prefix of each case: a block is evidence for
# one reference (CVS and the verification matrix take one per test).
lrt_result <- local({
  result <- NULL
  function() {
    if (is.null(result)) {
      source(file.path("lrt", "common.R"), local = TRUE)
      result <<- lrt_comparison("lrt")
    }
    result
  }
})

test_that("the likelihood-ratio kit compares every case", {
  skip_on_cran()
  result <- lrt_result()
  expect_identical(nrow(result), 27L)
  expect_true(all(grepl("^(mclust|glca|Mplus) ", result$case)), info = paste(unique(result$case), collapse = "; "))
})

invisible(lapply(list(list(name = "mclust and glca", pattern = "^(mclust|glca) "),
                      list(name = "Mplus", pattern = "^Mplus ")), function(program) {
  test_that(sprintf("likelihood-ratio tests agree with %s", program$name), {
    skip_on_cran()
    result <- lrt_result()
    rows <- result[grepl(program$pattern, result$case), , drop = FALSE]
    expect_gt(nrow(rows), 0L)
    disagreeing <- subset(rows, !agrees)
    expect_identical(nrow(disagreeing), 0L,
                     info = paste(disagreeing$case, disagreeing$quantity, collapse = "; "))
  })
}))
