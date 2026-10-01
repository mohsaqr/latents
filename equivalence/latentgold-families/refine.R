# Node refinement: if the base-run gaps are Latent GOLD's quadrature error,
# rerunning it with more nodes must move it towards latents. Compares each
# refined listing (returned/refine/) with latents and with the base run.
# Run from the project root.
devtools::load_all(".", quiet = TRUE)
base <- file.path("equivalence", "latentgold-families")
cases <- readRDS(file.path(base, "kit", "cases.rds"))
refined <- sub("\\.lst$", "", list.files(file.path(base, "returned", "refine"),
                                          pattern = "^g[0-9]+\\.lst$"))
read_lg <- function(directory, name, H) {
  lines <- iconv(readLines(file.path(directory, paste0(name, ".lst")), warn = FALSE),
                 "latin1", "UTF-8")
  ll <- as.numeric(strsplit(grep("Log-likelihood (LL)", lines, fixed = TRUE,
                                 value = TRUE)[1L], "\t")[[1L]][2L])
  posteriors <- utils::read.delim(file.path(directory, paste0(name, "_posteriors.txt")),
                                  check.names = FALSE, fileEncoding = "latin1")
  columns <- grep("^GClass#[0-9]+$", names(posteriors), value = TRUE)
  list(ll = ll, posterior = unique(cbind(g = posteriors[[1L]], posteriors[columns])))
}
gap <- function(mine, lg, H) {
  theirs <- as.matrix(lg[, -1L])
  mine <- mine[as.character(lg$g), , drop = FALSE]
  options <- if (H == 2L) list(1:2, 2:1) else
    list(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  min(vapply(options, function(o) max(abs(mine - theirs[, o])), numeric(1)))
}
rows <- lapply(refined, function(name) {
  case <- cases[[name]]
  data <- utils::read.delim(file.path(base, "kit", paste0(name, ".dat")))
  fit <- multilpa(data, paste0("y", seq_len(case$d)), "g", n_group_classes = case$H,
                  family = case$family, between_variance = case$between,
                  tol = 1e-14, max_iter = 20000, n_starts = 20, seed = 1)
  coarse <- read_lg(file.path(base, "returned"), name, case$H)
  fine <- read_lg(file.path(base, "returned", "refine"), name, case$H)
  data.frame(case = name, latents_ll = fit$log_likelihood,
             ll_gap_base = abs(fit$log_likelihood - coarse$ll),
             ll_gap_refined = abs(fit$log_likelihood - fine$ll),
             posterior_gap_base = gap(fit$group_posteriors, coarse$posterior, case$H),
             posterior_gap_refined = gap(fit$group_posteriors, fine$posterior, case$H))
})
result <- do.call(rbind, rows)
print(result, digits = 4, row.names = FALSE)
utils::write.csv(result, file.path(base, "refine.csv"), row.names = FALSE)
