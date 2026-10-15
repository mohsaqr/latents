test_that("pooled and bootstrap plot calls draw once and return reusable ggplots", {
  skip_if_not_installed("ggplot2")
  set.seed(192)
  d <- data.frame(g = rep(1:25, each = 4), a = rnorm(100), b = rnorm(100))
  pooled <- pool_imputations(list(d, d), c("a", "b"), "g", 1,
                             n_group_classes = 1, n_starts = 1, seed = 1)
  bootstrap <- structure(list(
    replicates = data.frame(statistic = seq(0, 8, length.out = 19)),
    statistic = 4, n_valid = 19L, iter = 19L, p_value = 0.55),
    class = "multilpa_bootstrap_lrt")
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)
  printed <- 0L
  plot_class <- class(ggplot2::ggplot())[1L]
  original <- utils::getS3method("print", plot_class)
  drawing <- function(x, ...) {
    printed <<- printed + 1L
    original(x, ...)
  }
  base::registerS3method("print", plot_class, drawing,
                         envir = asNamespace("ggplot2"))
  on.exit(base::registerS3method("print", plot_class, original,
                                 envir = asNamespace("ggplot2")), add = TRUE)
  invisible(lapply(seq_len(2L), function(i) plot(pooled, what = "measurement")))
  expect_equal(printed, 2L)
  returned <- withVisible(plot(bootstrap))
  expect_equal(printed, 3L)
  expect_false(returned$visible)
  expect_s3_class(returned$value, "ggplot")
  expect_s3_class(returned$value + ggplot2::theme(legend.position = "bottom"), "ggplot")
  expect_gt(length(grid::grid.ls(print = FALSE)$name), 0)
})
