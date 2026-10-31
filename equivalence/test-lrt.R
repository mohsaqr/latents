# compare_models(test = "lrt") against mclust (covariance structures), glca's
# gofglca() (a membership covariate) and Mplus (LTA measurement invariance).
# The references are stored in lrt/reference.rds by lrt/make-kit.R.
test_that("likelihood-ratio tests agree with mclust, glca and Mplus", {
  skip_on_cran()
  source(file.path("lrt", "common.R"), local = TRUE)
  result <- lrt_comparison("lrt")
  expect_identical(nrow(result), 27L)
  disagreeing <- subset(result, !agrees)
  expect_identical(nrow(disagreeing), 0L,
                   info = paste(disagreeing$case, disagreeing$quantity, collapse = "; "))
})
