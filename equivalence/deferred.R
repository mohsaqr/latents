# Load a deferred feature from `future/R/` for a validation script.
#
# `fit_random_intercept()` and `lmr_lrt()` were moved out of the shipped package
# (see `future/README.md`), but the comparisons that validated them are the
# evidence needed to bring them back, so they must keep running. The deferred
# file is evaluated in an environment whose parent is the package namespace:
# its functions see the current internal helpers. Its S3 methods are registered
# for the current session so generic calls can find methods enclosed here;
# no source, exports or installed package files are changed.
#
# Source this after `pkgload::load_all(".")`, then e.g.
#   lmr_lrt <- deferred_verb("lmr.R", "lmr_lrt")

#' One function from a deferred source file
#'
#' @param file Basename of a file in `future/R/`.
#' @param name The function to return from it.
#' @param root Project root directory.
#' @return The function, with its enclosure's parent set to the latents
#'   namespace.
deferred_verb <- function(file, name, root = ".") {
  stopifnot(
    "latents must be loaded first" = "latents" %in% loadedNamespaces(),
    "`file` must be a single file name" = is.character(file) && length(file) == 1L,
    "`name` must be a single function name" = is.character(name) && length(name) == 1L)
  path <- file.path(root, "future", "R", file)
  if (!file.exists(path)) {
    stop(errorCondition(sprintf("No deferred source at `%s`.", path),
                        class = "multilpa_missing_deferred", call = NULL))
  }
  environment <- new.env(parent = asNamespace("latents"))
  sys.source(path, envir = environment, keep.source = FALSE)
  generics <- c("as.data.frame", "parameter_inference", "get_results", "print",
                 "summary", "plot", "coef", "vcov", "confint", "logLik", "nobs")
  invisible(lapply(ls(environment, all.names = TRUE), function(method) {
    generic <- generics[startsWith(method, paste0(generics, "."))]
    if (length(generic) == 1L && is.function(environment[[method]])) {
      registerS3method(generic, substring(method, nchar(generic) + 2L),
                       environment[[method]], envir = asNamespace("latents"))
    }
  }))
  if (!exists(name, envir = environment, inherits = FALSE)) {
    stop(errorCondition(sprintf("`%s` does not define `%s()`.", path, name),
                        class = "multilpa_missing_deferred", call = NULL))
  }
  get(name, envir = environment, inherits = FALSE)
}
