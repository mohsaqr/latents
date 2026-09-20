# Every test that draws must draw onto a device it opened itself. Base graphics
# always has an active device; when no test opened one, R silently opens
# `pdf("Rplots.pdf")` in the working directory, so a forgotten device does not
# fail the test -- it leaves a file in tests/testthat instead.
draw <- function(expression) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 900, height = 600)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE, after = FALSE)
  force(expression)
}
