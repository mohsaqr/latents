# The three-step latent Markov model with group classes (r3step() on
# lta(n_group_classes = 2)) against Latent GOLD 6.1's `step3 modal ml` with a
# static GClass beside the dynamic State. PENDING: the kit in
# latentgold-step3-markov-groups/kit was built on a machine without Latent
# GOLD, so returned/ holds no listings yet and this test skips. When the
# listings are returned, the log likelihood is compared here at Latent GOLD's
# printed precision; the coefficient table's layout with a second
# conditioning latent variable is to be read from the first listing and its
# parser added (see the kit's README).

base <- file.path("latentgold-step3-markov-groups")
cases <- readRDS(file.path(base, "kit", "cases.rds"))

invisible(lapply(names(cases), function(name) {
  test_that(sprintf("group-class three-step latent Markov model matches Latent GOLD: %s",
                    name), {
    listing <- file.path(base, "returned", paste0(name, ".lst"))
    skip_if_not(file.exists(listing),
                "Latent GOLD listing pending (no Latent GOLD on the build machine)")
    case <- cases[[name]]
    fit <- readRDS(file.path(base, "kit", paste0(name, "-fit.rds")))
    data <- readRDS(file.path(base, "kit", paste0(name, "-data.rds")))
    ours <- r3step(fit, data, initial_covariates = case$initial_covariates,
                   transition_covariates = case$transition_covariates,
                   vcov_type = "observed")
    lines <- sub("\r$", "", iconv(readLines(listing, warn = FALSE), "latin1", "UTF-8"))
    row <- grep("^Log-likelihood \\(LL\\)", lines, value = TRUE)[1L]
    expect_lt(abs(attr(ours, "log_likelihood") -
                    as.numeric(strsplit(row, "\t")[[1L]][2L])), 5.1e-5)
  })
}))
