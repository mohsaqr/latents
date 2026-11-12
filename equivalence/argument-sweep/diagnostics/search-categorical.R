# How often does each search reach the best known maximum of a latent class
# model? Ten starts of the default search (k-means), of emEM, and of poLCA's
# random starts, against the best of every run including 40-start searches.
# (With 100-start arms and five seeds the study ran past ninety minutes on a
# shared machine and was stopped.)
Sys.setenv(NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("equivalence/argument-sweep/helpers.R")
source("equivalence/argument-sweep/datasets.R")
environment <- new.env()
utils::data(list = c("cheating", "gss82", "election"), package = "poLCA", envir = environment)
as_factors <- function(data) as.data.frame(lapply(data, factor))
datasets <- c(sweep_categorical(), list(
  cheating = as_factors(environment$cheating[1:4]), gss82 = as_factors(environment$gss82),
  election = as_factors(environment$election[1:12]),
  esm_missing = sweep_punch(sweep_categorical()$esm, 0.1)))
cases <- sweep_grid(dataset = names(datasets), k = 2:3, seed = 1:3)
rows <- parallel::mclapply(cases, function(o) {
  data <- datasets[[o$dataset]]
  codes <- as.data.frame(lapply(data, function(v) as.integer(factor(v))))
  formula <- stats::as.formula(sprintf("cbind(%s) ~ 1", paste(names(data), collapse = ", ")))
  ours <- function(search, n_starts) {
    suppressWarnings(lca(data, names(data), n_classes = o$k, missing = "fiml", seed = o$seed,
                         n_starts = n_starts, search = search, tol = 1e-10,
                         max_iter = 20000L))$log_likelihood
  }
  polca <- function(nrep) sweep_seeded(o$seed, function() {
    suppressWarnings(poLCA::poLCA(formula, codes, nclass = o$k, maxiter = 20000L, tol = 1e-10,
                                  na.rm = FALSE, verbose = FALSE, calc.se = FALSE,
                                  nrep = nrep))$llik
  })
  found <- c(kmeans_10 = ours("kmeans", 10L), emem_10 = ours("emEM", 10L), polca_10 = polca(10L),
             emem_40 = ours("emEM", 40L), polca_40 = polca(40L))
  data.frame(dataset = o$dataset, k = o$k, seed = o$seed, best = max(found),
             t(max(found) - found))
}, mc.cores = 6L)
table <- do.call(rbind, rows)
options(width = 200, digits = 6)
reached <- function(gap) mean(gap < 1e-3)
summary <- stats::aggregate(cbind(kmeans_10, emem_10, polca_10, emem_40, polca_40) ~ dataset + k,
                            table, reached)
cat("Share of three seeds reaching the best known log likelihood (within 1e-3):\n")
print(summary, row.names = FALSE)
cat("\nOverall share reaching the best known maximum:\n")
print(round(vapply(table[c("kmeans_10", "emem_10", "polca_10", "emem_40", "polca_40")],
                   reached, numeric(1)), 3))
cat("\nMean shortfall in log likelihood:\n")
print(round(colMeans(table[c("kmeans_10", "emem_10", "polca_10", "emem_40", "polca_40")]), 4))
