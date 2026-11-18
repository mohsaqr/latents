# Exports Latent GOLD comparison cases as bundles for Prova Verify (format:
# plan/24-verify.md in the prova repository), with Prova's R writer. A bundle
# sets latents beside Latent GOLD: the code of each, the raw output of each,
# the data, and every compared quantity with a value per side.
#
# The latents side is an R script derived from the kit's own fitting call. It
# is stored in the bundle and run there in its own R process, so the code and
# the console shown beside it are one run. The Latent GOLD side is the syntax
# of the kit and the output retained from the licensed run. The compared
# numbers are the committed `comparison.csv` rows; a bundle is refused when the
# refitted log likelihood is not the one recorded there.
#
# Run from the project root:
#   Rscript equivalence/latentgold/export-bundle.R [case ...]
devtools::load_all(".", quiet = TRUE)
source(Sys.getenv("PROVA_VERIFY_R", file.path("..", "prova", "reporters", "verify.R")))
source(file.path("equivalence", "latentgold", "cases.R"))

.bundle_base <- file.path("equivalence", "latentgold")
.bundle_generator <- "latents equivalence/latentgold/export-bundle.R"

#' The model call of a case, as it would be typed
#'
#' The kit holds each case's fit as a closure, `.lg_fit(<data>, ...)` or a block
#' whose first statement is `fit <- lta(<data>, ...)`. The data argument becomes
#' `data`, every other bare name is replaced by its value in the closure, and
#' `.lg_fit()` is written out as the `multilpa()` call it makes.
#' @param case One element of [lg_cases()].
#' @return The call, or `NULL` when the case's fit is not one of the two forms.
.bundle_call <- function(case) {
  body <- body(case$fit)
  values <- environment(case$fit)
  call <- if (is.call(body) && identical(body[[1L]], as.name("{"))) {
    first <- body[[2L]]
    if (is.call(first) && identical(first[[1L]], as.name("<-"))) first[[3L]]
  } else body
  if (!is.call(call) || !as.character(call[[1L]]) %in% c(".lg_fit", "lta")) return(NULL)
  arguments <- as.list(call)[-1L]
  arguments[[1L]] <- as.name("data")
  arguments[-1L] <- lapply(arguments[-1L], function(argument) {
    if (is.name(argument)) get(as.character(argument), envir = values) else argument
  })
  if (identical(call[[1L]], as.name(".lg_fit"))) {
    arguments <- c(arguments, list(n_starts = 50L, max_iter = 200000L, tol = 1e-14, seed = 1L))
    as.call(c(as.name("multilpa"), arguments))
  } else {
    as.call(c(as.name("lta"), arguments))
  }
}

#' The latents script of a case
#'
#' @param case One element of [lg_cases()].
#' @return The script, one line per element, or `NULL` when the case's model
#'   call cannot be derived ([.bundle_call()]).
.bundle_script <- function(case) {
  call <- .bundle_call(case)
  if (is.null(call)) return(NULL)
  fitted <- deparse(call, width.cutoff = 80L)
  c("library(latents)",
    # The kit writes a missing value as ".", Latent GOLD's code for it.
    sprintf("data <- utils::read.delim(\"../data/%s\", na.strings = \".\")", case$data_file),
    paste0(c("fit <- ", rep("       ", length(fitted) - 1L)), trimws(fitted)),
    "print(fit)",
    "summary(fit)",
    "result <- data.frame(quantity = \"log likelihood\", value = as.numeric(logLik(fit)))")
}

#' The function a script fits with, as the case's `covers`
#'
#' @param script A script made by [.bundle_script()].
#' @return `"multilpa()"` or `"lta()"`.
.bundle_covers <- function(script) {
  fitted <- grep("^fit <- ", script, value = TRUE)
  sprintf("%s()", sub("^fit <- ([A-Za-z_.]+)\\(.*$", "\\1", fitted))
}

#' What a compared quantity is, from the name the comparison gives it
#'
#' @param quantity The quantity names of `comparison.csv`.
#' @return One word or phrase per quantity.
.bundle_kind <- function(quantity) {
  kinds <- c("^log likelihood" = "likelihood", "^parameters$|^n_parameters" = "parameter count",
             "^n_|^number of" = "sample size", "^(aic|bic|caic|sabic)" = "information criterion",
             "^posteriors" = "posterior", "^standard error" = "standard error",
             "^(mean|variance|estimate)|proportion$" = "estimate",
             "residual" = "residual")
  hits <- vapply(names(kinds), function(pattern) grepl(pattern, quantity, ignore.case = TRUE),
                 logical(length(quantity)))
  first <- apply(matrix(hits, nrow = length(quantity)), 1L, function(row) match(TRUE, row))
  ifelse(is.na(first), "", unname(kinds)[first])
}

#' The compared quantities of one case, as manifest rows
#'
#' @param name The case's name.
#' @param listing,posteriors The bundle paths of the Latent GOLD listing and
#'   classification file, named as the source of a reference value by the
#'   comparison's own note.
#' @return A list with one element per compared quantity: the quantity, what
#'   it is, the value of each side, the difference, the tolerance and the status.
.bundle_comparisons <- function(name, listing, posteriors) {
  rows <- subset(utils::read.csv(file.path(.bundle_base, "comparison.csv")), case == name)
  stopifnot("the case has no rows in comparison.csv" = nrow(rows) > 0L)
  source <- ifelse(grepl("posteriors", rows$note, fixed = TRUE), posteriors, listing)
  Map(function(quantity, kind, reference, obtained, difference, tolerance, status, note, path) {
    list(quantity = quantity, kind = kind,
         values = list(latentgold = reference, latents = obtained),
         baseline = "latentgold", difference = difference, measure = "absolute",
         tolerance = tolerance, status = status, note = note,
         sources = list(latentgold = list(path = path)))
  }, rows$quantity, .bundle_kind(rows$quantity), rows$reference, rows$obtained, rows$difference,
  rows$tolerance, gsub(" ", "_", rows$status, fixed = TRUE), rows$note, source) |> unname()
}

#' Export one case
#'
#' @param name The case's name, as in `comparison.csv`.
#' @param out The folder that receives the bundle folder.
#' @return One row: the case, whether it was exported, how many quantities it
#'   compares, and the reason when it was not exported.
export_bundle <- function(name, out = file.path(.bundle_base, "bundles")) {
  cases <- lg_cases()
  stopifnot("`name` must be one case of the kit" = length(name) == 1L && name %in% names(cases))
  case <- cases[[name]]
  skipped <- function(reason) data.frame(case = name, exported = FALSE, compared = NA_integer_, reason = reason)
  script <- .bundle_script(case)
  if (is.null(script)) return(skipped("the model call is not a single multilpa() or lta() call"))

  bundle <- file.path(out, name)
  unlink(bundle, recursive = TRUE)
  dir.create(bundle, recursive = TRUE)
  kit <- file.path(.bundle_base, "kit")
  returned <- file.path(.bundle_base, "returned")
  latents_side <- verify_side("latents", package = "latents", primary = TRUE, load = ".", script = script)
  latentgold_side <- verify_side(
    "latentgold", name = "Latent GOLD", version = "6.1", authority = "commercial", source = "retained",
    run_at = "2026-10-19", platform = "Wine Devel 11.17 on macOS 26.3.1, arm64",
    licence = "free academic single-user licence", command = "lg61.exe <case>.lgs /b /o <case>.lst")

  data <- verify_copy(file.path(kit, case$data_file), bundle, file.path("data", case$data_file),
                      label = "Data", type = "table", rows = nrow(case$data), columns = names(case$data),
                      origin = "From the comparison kit, equivalence/latentgold/cases.R.")
  from_kit <- function(files, folder, ...) {
    lapply(files[file.exists(file.path(folder, files))], function(file)
      verify_copy(file.path(folder, file), bundle, file.path("latentgold", file), side = "latentgold", ...))
  }
  models <- sub("\\.lgs$", "", names(case$files))
  syntax <- from_kit(paste0(models, ".lgs"), kit, label = "Syntax", type = "code", language = "latentgold")
  listings <- from_kit(paste0(models, ".lst"), returned, label = "Listing", type = "listing", encoding = "latin1")
  classification <- from_kit(paste0(models, "_posteriors.txt"), returned, label = "Classification file",
                             type = "table", encoding = "latin1")
  if (length(listings) == 0L) return(skipped("no Latent GOLD listing is retained for it"))

  run <- verify_run(bundle, latents_side, reads = data$path)
  comparisons <- .bundle_comparisons(name, listing = listings[[1L]]$path,
                                     posteriors = c(classification, listings)[[1L]]$path)
  recorded <- Filter(function(row) identical(row$quantity, "log likelihood"), comparisons)[[1L]]$values$latents
  if (!isTRUE(all.equal(run$result$value, recorded, tolerance = 1e-10))) {
    stop(errorCondition(sprintf(
      "Case `%s`: the refitted log likelihood %.10f is not the %.10f in comparison.csv; rerun compare.R first.",
      name, run$result$value, recorded), class = "latents_stale_bundle_numbers", call = NULL))
  }
  sidecar <- subset(utils::read.csv(file.path(.bundle_base, "comparison.csv.provenance.csv")), kind == "result")

  verify_write(
    bundle, title = name, description = case$description, covers = .bundle_covers(script),
    sides = list(latents_side, latentgold_side),
    sections = list(
      list(id = "code", title = "Code", items = c(list(run$code), syntax)),
      list(id = "output", title = "Raw output", items = c(list(run$output), listings, classification)),
      list(id = "data", title = "Data", items = list(data))),
    comparisons = comparisons, generator = .bundle_generator, id = paste0("latentgold/", name),
    comparison_source = list(file = "equivalence/latentgold/comparison.csv",
                             version = sidecar$version, recorded_at = sidecar$recorded_at))
  data.frame(case = name, exported = TRUE, compared = length(comparisons), reason = "")
}

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- names(lg_cases())
print(do.call(rbind, lapply(requested, export_bundle)), row.names = FALSE)
