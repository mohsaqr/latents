# Latent GOLD 6.1 against the package as it is now (latentgold/).
#
# latentgold/compare.R compares Latent GOLD's retained output with targets the
# package wrote when the kit was built, so it cannot catch a later regression
# in the package. This test builds the targets from the current package,
# case by case, and applies the same comparison: likelihood, parameters,
# criteria, posteriors at each level, means and variances implied by Latent
# GOLD's posteriors, standard errors where both report them, and the
# three-step estimates. Any quantity that disagrees, or whose Latent GOLD
# output is missing or unreadable, fails the test.

latentgold_root <- file.path("latentgold")
local({
  source(file.path(latentgold_root, "cases.R"), local = TRUE)
  source(file.path(latentgold_root, "compare-functions.R"), local = TRUE)
  source(file.path(latentgold_root, "targets.R"), local = TRUE)

  test_that("every Latent GOLD case agrees with the current package", {
    cases <- lg_cases()
    comparison <- do.call(rbind, lapply(names(cases), function(name) {
      case <- cases[[name]]
      target <- .targets(name, case, suppressWarnings(case$fit()))
      lg_compare_case(target, file.path(latentgold_root, "returned"))
    }))
    disagreeing <- subset(comparison, status == "disagree")
    expect_identical(nrow(disagreeing), 0L,
                     info = paste(capture.output(print(disagreeing)), collapse = "\n"))
    absent <- subset(comparison, status %in% c("missing", "unparsed"))
    expect_identical(nrow(absent), 0L,
                     info = paste(capture.output(print(absent)), collapse = "\n"))
    expect_gt(sum(comparison$status == "agree"), 300L)
  })
})
