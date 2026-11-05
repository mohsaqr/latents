# Mutation check of the area's central formulas. Applies one edit at a time to a
# throwaway copy of the package (never the repo), and runs the single test files
# that name that formula. Usage: Rscript this.R <first> <last>.
args <- as.numeric(commandArgs(TRUE)[1:2])
extra_files <- commandArgs(TRUE)[-(1:2)]
repo <- "/Users/mohammedsaqr/Documents/Github/latents"
base <- "/private/tmp/claude-503/-Users-mohammedsaqr-Documents-Github-latents/2f29606d-312f-4450-b6bf-33795875bed9/scratchpad/mut"
M <- list(
 list("KIC 3(q+1) -> 3q", "diagnostics.R", "penalty = c(0, 2 * q, 3 * (q + 1))", "penalty = c(0, 2 * q, 3 * q)", c("test-diagnostics.R")),
 list("SABIC (n+2)/24 -> n/24", "diagnostics.R", "sabic = q * log((n + 2) / 24)", "sabic = q * log(n / 24)", c("test-diagnostics.R")),
 list("CLC 2EN -> EN", "diagnostics.R", "clc = 2 * entropy_sum)", "clc = entropy_sum)", c("test-diagnostics.R","test-enumeration.R")),
 list("AWE 1.5 -> 1", "diagnostics.R", "awe = 2 * q * (1.5 + log(n)) + 2 * entropy_sum", "awe = 2 * q * (1 + log(n)) + 2 * entropy_sum", c("test-diagnostics.R")),
 list("CAIC +1 dropped", "diagnostics.R", "caic = q * (log(n) + 1)", "caic = q * log(n)", c("test-diagnostics.R")),
 list("OCC proportion odds dropped", "diagnostics.R", "(proportion[defined] / (1 - proportion[defined]))", "proportion[defined]", c("test-diagnostics.R","test-review-1027-diagnostics.R")),
 list("BCH rows<->columns of D^-1", "three-step.R", "weights <- inverse[pieces$modal, , drop = FALSE]", "weights <- t(inverse)[pieces$modal, , drop = FALSE]", c("test-three-step.R")),
 list("error matrix normalised by assigned", "three-step.R", "colSums(weighted * (pieces$modal == assigned)) / totals", "colSums(weighted * (pieces$modal == assigned)) / sum(weighted * (pieces$modal == assigned))", c("test-three-step.R")),
 list("aLMR drops d", "lmr.R", "adjusted <- lr / (1 + 1 / (added * log(n_units)))", "adjusted <- lr / (1 + 1 / log(n_units))", c("test-lmr.R","test-features-vlmr.R")),
 list("Vuong null block sign", "lmr.R", "signs <- c(rep(1, ncol(alternative$scores)), rep(-1, ncol(null$scores)))", "signs <- c(rep(1, ncol(alternative$scores)), rep(1, ncol(null$scores)))", c("test-lmr.R","test-features-vlmr.R")),
 list("BLRT p (1+b)/(V+1) -> b/V", "simulation-service.R", "(1 + sum(replicates$statistic[valid] >= observed)) / (n_valid + 1)", "sum(replicates$statistic[valid] >= observed) / n_valid", c("test-bootstrap-lrt-extended.R","test-enumeration.R")),
 list("knee distance not perpendicular", "selection.R", "distance <- below / sqrt(1 + slope^2)", "distance <- below", c("test-enumeration.R")),
 list("Schwarz exp(-d/2) -> exp(-d)", "compare-models.R", "relative <- exp(-table$delta_bic / 2)", "relative <- exp(-table$delta_bic)", c("test-enumeration.R","test-compare-lrt.R")),
 list("BVR Fisher n-3 -> n", "bivariate-residuals.R", "* sqrt(effective - 3)", "* sqrt(effective)", c("test-bivariate-residuals.R")),
 list("BVR df RC-1", "bivariate-residuals.R", "degrees <- (nrow(observed) - 1L) * (ncol(observed) - 1L)", "degrees <- nrow(observed) * ncol(observed) - 1L", c("test-bivariate-residuals.R")),
 list("tree flow unweighted", "tree.R", "crossprod(posteriors[[row - 1L]] * weights, posteriors[[row]]) / total", "crossprod(posteriors[[row - 1L]], posteriors[[row]]) / total", c("test-tree.R")),
 list("relative entropy log(K+1)", "diagnostics.R", "(total * log(ncol(probabilities)))", "(total * log(ncol(probabilities) + 1))", c("test-diagnostics.R")),
 list("EN drops group-class entropy", "classification-entropy.R", "row_part + .multilpa_entropy_sum(", "row_part + 0 * .multilpa_entropy_sum(", c("test-diagnostics.R")),
 list("simulation ignores group class", "enumeration.R", "object$profile_probabilities[group_class[object$group_index], , drop = FALSE])", "object$profile_probabilities[rep(1L, object$n_observations), , drop = FALSE])", c("test-bootstrap-lrt-extended.R","test-enumeration.R")),
 list("BVR expected uses wrong profile total", "bivariate-residuals.R", "replace(numeric(object$n_profiles), profile, sum(pair_weight))", "replace(numeric(object$n_profiles), profile, sum(weight))", c("test-bivariate-residuals.R")),
 list("LRT chi-square df+1", "compare-lrt.R", "df <- as.integer(alternative$n_parameters - null$n_parameters)", "df <- as.integer(alternative$n_parameters - null$n_parameters + 1L)", c("test-compare-lrt.R")),
 list("LMR step from one more profile also nested group", "alignment-service.R", "moved - moved[, if (identical(reference, \"last\")) ncol(moved) else 1L]", "moved", c("test-alignment-service.R","test-align-labels.R"))
)
run <- function(i) {
  m <- M[[i]]; dir <- file.path(base, sprintf("m%02d", i)); unlink(dir, recursive = TRUE); dir.create(dir)
  invisible(file.copy(file.path(repo, c("DESCRIPTION","NAMESPACE","R","data","inst","tests")), dir, recursive = TRUE))
  path <- file.path(dir, "R", m[[2]]); src <- readLines(path); joined <- paste(src, collapse = "\n")
  stopifnot("mutation text not found" = grepl(m[[3]], joined, fixed = TRUE))
  writeLines(strsplit(gsub(m[[3]], m[[4]], joined, fixed = TRUE), "\n", fixed = TRUE)[[1]], path)
  suppressMessages(devtools::load_all(dir, quiet = TRUE, helpers = FALSE))
  files <- unique(c(m[[5]], extra_files))
  out <- lapply(files, function(f) {
    fp <- file.path(dir, "tests", "testthat", f)
    if (!file.exists(fp)) return(data.frame(file = f, failed = NA, passed = NA))
    r <- as.data.frame(suppressMessages(testthat::test_file(fp, reporter = "silent", stop_on_failure = FALSE, env = new.env(parent = asNamespace("latents")))))
    data.frame(file = f, failed = sum(r$failed) + sum(r$error), passed = sum(r$passed))
  })
  o <- do.call(rbind, out); unlink(dir, recursive = TRUE)
  cat(sprintf("%-46s %s\n", m[[1]], paste(sprintf("%s: %s fail/%s pass", o$file, o$failed, o$passed), collapse = "; ")))
}
Sys.setenv(NOT_CRAN = "true"); pdf(NULL); set.seed(1)
invisible(lapply(seq(args[1], args[2]), run))
