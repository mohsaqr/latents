# Writes test-datasets.csv: one row per dataset that a testthat file in
# equivalence/ compares with another program, read from the file's own
# simulation kit. Only datasets the external program actually ran are listed
# (the VLMR kit's 62 of 107, the ML three-step kit's 70 of 193); datasets
# checked against independent code alone are not. Run from the project root
# whenever a kit gains or loses cases:
#   Rscript equivalence/verification-matrix/test-datasets.R

#' One row per dataset: the test file, the function it checks, the dataset
#' @noRd
.dataset_rows <- function(file, fn, kit, names) {
  if (!length(names)) return(NULL)
  data.frame(location = file.path("equivalence", file), fn = fn,
             dataset = paste0(kit, ": ", names), stringsAsFactors = FALSE)
}

#' The cases of a kit, named (by their `name`, or their list names, or their
#' position)
#' @noRd
.kit_cases <- function(path) {
  cases <- readRDS(file.path("equivalence", path))
  labels <- vapply(seq_along(cases), function(i) {
    cases[[i]]$name %||% names(cases)[i] %||% sprintf("case %03d", i)
  }, character(1))
  labels[!nzchar(labels)] <- sprintf("case %03d", which(!nzchar(labels)))
  stats::setNames(cases, labels)
}

#' Every row of test-datasets.csv
#' @noRd
test_datasets <- function() {
  montecarlo <- .kit_cases(file.path("montecarlo", "kit", "cases.rds"))
  montecarlo_fn <- vapply(montecarlo, function(case) {
    categorical <- !nzchar(case$spec$model %||% "")
    if (identical(case$spec$level, "two")) {
      if (categorical) "multilca" else "multilpa"
    } else if (categorical) "lca" else "lpa"
  }, character(1))
  extended <- .kit_cases(file.path("montecarlo", "extended", "kit", "cases.rds"))
  designs <- readRDS(file.path("equivalence", "montecarlo", "enumeration", "designs.rds"))
  vlmr <- .kit_cases(file.path("features", "vlmr", "kit", "cases.rds"))
  ml3 <- .kit_cases(file.path("features", "ml-three-step", "kit", "cases.rds"))
  pair <- .kit_cases(file.path("features", "pair-covariance", "kit", "cases.rds"))
  joint <- .kit_cases(file.path("features", "joint-items", "kit", "cases.rds"))
  joint_fn <- ifelse(vapply(joint, function(case) is.null(case$spec$h), logical(1)),
                      "lca", "multilca")
  direct <- .kit_cases(file.path("features", "direct-effects", "kit", "cases.rds"))
  families <- .kit_cases(file.path("latentgold-families", "kit", "cases.rds"))
  markov <- .kit_cases(file.path("latentgold-step3-markov", "kit", "cases.rds"))
  missing <- .kit_cases(file.path("missing-patterns", "kit", "cases.rds"))
  missing_fn <- ifelse(vapply(missing, function(case) identical(case$level, "two"), logical(1)),
                       "multilpa", "lpa")
  external_vlmr <- vapply(vlmr, function(case) isTRUE(case$external), logical(1))
  external_ml3 <- vapply(ml3, function(case) isTRUE(case$lg) || isTRUE(case$mplus), logical(1))
  external_pair <- vapply(pair, function(case) length(case$spec$programs) > 0L, logical(1))
  rows <- c(
    lapply(unique(montecarlo_fn), function(fn) {
      .dataset_rows("test-montecarlo.R", fn, "Monte Carlo kit", names(montecarlo)[montecarlo_fn == fn])
    }),
    list(
      .dataset_rows("test-montecarlo-extended.R", "multilpa", "Monte Carlo extended kit", names(extended)),
      .dataset_rows("test-enumeration-mplus.R", "enumeration", "Monte Carlo kit", names(designs)),
      .dataset_rows("test-features-vlmr.R", "bootstrap_lrt", "VLMR kit", names(vlmr)[external_vlmr]),
      .dataset_rows("test-features-ml-three-step.R", "three_step", "ML three-step kit",
                    names(ml3)[external_ml3]),
      .dataset_rows("test-features-pair-covariance.R", "lpa", "pair-covariance kit",
                    names(pair)[external_pair]),
      .dataset_rows("test-features-direct-effects.R", "lpa", "direct-effects kit", names(direct)),
      .dataset_rows("test-latentgold-families.R", "multilpa", "Latent GOLD families kit", names(families)),
      .dataset_rows("test-latentgold-step3-markov.R", "three_step", "Latent GOLD step-3 Markov kit",
                    names(markov))),
    lapply(unique(joint_fn), function(fn) {
      .dataset_rows("test-features-joint-items.R", fn, "joint-items kit", names(joint)[joint_fn == fn])
    }),
    lapply(unique(missing_fn), function(fn) {
      .dataset_rows("test-missing-patterns.R", fn, "missing-patterns kit", names(missing)[missing_fn == fn])
    }))
  out <- do.call(rbind, rows)
  stopifnot("no dataset is listed twice for a file and function" =
              !anyDuplicated(out[c("location", "fn", "dataset")]))
  out
}

datasets <- test_datasets()
utils::write.csv(datasets, file.path("equivalence", "verification-matrix", "test-datasets.csv"),
                 row.names = FALSE)
print(stats::aggregate(dataset ~ location + fn, datasets, length), row.names = FALSE)
