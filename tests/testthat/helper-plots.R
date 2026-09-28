# A plot() method returns a ggplot object rather than drawing. Building it
# runs every layer's computation, which is where a broken view fails, so these
# expectations build the plot and fail on any warning.
expect_plot <- function(object) {
  testthat::expect_s3_class(object, "ggplot")
  testthat::expect_no_warning(ggplot2::ggplot_build(object))
  invisible(object)
}

# A list of plots, as from what = "all" or plot() on a diagnostics result.
expect_plots <- function(object) {
  testthat::expect_s3_class(object, "latents_plots")
  testthat::expect_gt(length(object), 0L)
  invisible(lapply(object, expect_plot))
  invisible(object)
}

# The data one layer of a built plot draws, for checking the numbers.
plot_layer <- function(object, layer) ggplot2::layer_data(object, layer)
