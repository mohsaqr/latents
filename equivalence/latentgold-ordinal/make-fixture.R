# Store the Latent GOLD results as a test fixture, so the test suite checks
# them without Wine. Run from the package root after compare.R's inputs exist.
base <- file.path("equivalence", "latentgold-ordinal")
cases <- readRDS(file.path(base, "kit", "cases.rds"))
fixture <- lapply(names(cases), function(name) {
  case <- cases[[name]]
  lines <- iconv(readLines(file.path(base, "returned", paste0(name, ".lst")), warn = FALSE),
                 "latin1", "UTF-8")
  grab <- function(label) {
    as.numeric(strsplit(grep(label, lines, value = TRUE)[1L], "\t")[[1L]][2L])
  }
  list(data = utils::read.delim(file.path(base, "kit", paste0(name, ".dat"))),
       vars = names(case$indicators),
       types = vapply(case$indicators, `[[`, character(1), "type"),
       n_profiles = length(case$sizes),
       n_group_classes = if (is.null(case$group_size)) 1L else length(case$group_classes),
       two_level = !is.null(case$group_size),
       latent_gold_ll = grab("^Log-likelihood \\(LL\\)"),
       latent_gold_npar = grab("^Number of parameters"))
})
names(fixture) <- names(cases)
saveRDS(fixture, file.path("tests", "testthat", "fixtures", "latentgold-ordinal.rds"),
        compress = "xz")
str(lapply(fixture, `[`, c("latent_gold_ll", "latent_gold_npar")))
