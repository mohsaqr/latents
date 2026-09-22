# Build the Latent GOLD comparison kit.
#
# Run from the project root:
#   Rscript validation/latentgold/make-kit.R
#
# Writes validation/latentgold/kit/: one tab-separated data file per dataset,
# one Latent GOLD syntax file per model, run-all.bat to estimate them in order,
# and MANIFEST.txt recording the multilpa version and an md5 for every file.
# The targets this package produced are written to
# validation/latentgold/targets/, where compare.R reads them.
#
# Copy kit/ to a Windows machine with Latent GOLD 6.1, run run-all.bat, and
# copy the .lst listings and *_posteriors.txt files it writes back into
# validation/latentgold/returned/. Then run compare.R.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("validation", "latentgold", "cases.R"))

directory <- file.path("validation", "latentgold")
kit <- file.path(directory, "kit")
targets_dir <- file.path(directory, "targets")
# The kit is regenerated whole, so a file from a case that no longer exists
# cannot survive into the next hand-off.
unlink(c(kit, targets_dir), recursive = TRUE)
dir.create(kit, recursive = TRUE)
dir.create(targets_dir, recursive = TRUE)

#' Posteriors of one level as one row per unit and one column per class
#'
#' @param long A long posterior table from `get_results()`.
#' @param unit The column naming the unit, `"row"` or `"group"`.
#' @param class The column naming the class, `"profile"` or `"group_class"`.
#' @return A data frame with `unit` and `class_1`, `class_2`, ...
.wide_posteriors <- function(long, unit, class) {
  classes <- sort(unique(long[[class]]))
  units <- unique(long[[unit]])
  columns <- lapply(classes, function(k) {
    rows <- long[[class]] == k
    long$posterior[rows][match(units, long[[unit]][rows])]
  })
  stopifnot("posterior table is not complete" = !anyNA(unlist(columns)))
  data.frame(unit = units, stats::setNames(columns, sprintf("class_%d", classes)))
}

#' Everything compare.R needs from this package for one case
#'
#' @param name The case's name.
#' @param case One element of [lg_cases()].
#' @param fit This package's fit to it.
#' @return A list of targets.
.targets <- function(name, case, fit) {
  # A target must be the optimum the search keeps finding, not one start's
  # luck; otherwise a disagreement with Latent GOLD is uninterpretable.
  starts <- get_results(fit, "starts")
  replicated <- sum(abs(starts$log_likelihood - fit$log_likelihood) < 1e-6,
                    na.rm = TRUE)
  if (replicated < 2L) {
    stop(errorCondition(sprintf(
      "Case `%s`: the best log likelihood was reached by %d start(s); a target needs at least 2.",
      name, replicated), class = "multilpa_unreplicated_target", call = NULL))
  }
  individual <- get_results(fit, "posteriors", format = "wide")
  profile_columns <- grep("^posterior_profile_", names(individual), value = TRUE)
  posteriors <- data.frame(id = case$data$id[individual$row],
                           stats::setNames(individual[profile_columns],
                                           sub("^posterior_profile_", "class_", profile_columns)))
  groups <- if (fit$n_group_classes > 1L) {
    wide <- .wide_posteriors(get_results(fit, "group_posteriors"), "group", "group_class")
    names(wide)[names(wide) == "unit"] <- "g"
    wide
  }
  inference <- if (name %in% c("c01_continuous_varying", "c02_continuous_equal",
                               "c07_covariates")) {
    parameter_inference(fit)
  }
  steps <- if (identical(name, "c09_three_step")) {
    list(classification_errors = get_results(fit, "classification_errors",
                                             level = "individuals"),
         distal_bch = three_step(fit, case$data, "distal", method = "bch"),
         distal_modal = three_step(fit, case$data, "distal", method = "modal"),
         distal_proportional = three_step(fit, case$data, "distal",
                                          method = "proportional"),
         covariate_ml = r3step(fit, case$data, "x"))
  }
  list(case = name, description = case$description, anchors = case$anchors,
       lgs = names(case$files), data_file = case$data_file,
       continuous = fit$continuous, variance_model = fit$variance_model,
       covariance_model = fit$covariance_model %||% "diagonal",
       missing = fit$missing %||% "error",
       log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters,
       n_best_replicated = replicated,
       posteriors = posteriors, group_posteriors = groups,
       profiles = get_results(fit, "profiles"),
       group_class_probabilities = if (fit$n_group_classes > 1L) fit$group_probabilities,
       parameters = inference,
       residuals = if (identical(name, "c08_bivariate_residuals"))
         get_results(fit, "residuals"),
       steps = steps,
       data = case$data,
       multilpa_version = as.character(utils::packageVersion("multilpa")))
}

#' Write lines with Windows line endings
#'
#' Everything in the kit is read on Windows, by Latent GOLD and by cmd.exe,
#' and a batch file with bare LF endings can misparse; so every kit file is
#' written with CRLF, not only the batch file.
#'
#' @param lines Character vector.
#' @param path File to write.
#' @return `path`, invisibly.
.write_windows <- function(lines, path) {
  # Evaluate the lines before the file exists: opening the connection creates
  # it, and a manifest computed lazily afterwards would list itself.
  force(lines)
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE, after = FALSE)
  writeLines(lines, connection, sep = "\r\n")
  invisible(path)
}

#' Write each dataset once, even when several cases read it
#'
#' Latent GOLD reads a header row and a tab-separated body; a period marks a
#' missing value.
#'
#' @param cases From [lg_cases()].
#' @param kit The kit directory.
#' @return The data file names, invisibly.
.write_kit_data <- function(cases, kit) {
  files <- vapply(cases, function(case) case$data_file, character(1))
  owners <- cases[!duplicated(files)]
  invisible(vapply(owners, function(case) {
    utils::write.table(case$data, file.path(kit, case$data_file), sep = "\t",
                       quote = FALSE, row.names = FALSE, na = ".", eol = "\r\n")
    case$data_file
  }, character(1)))
}

#' Every syntax file of every case, in estimation order
#' @param cases From [lg_cases()].
#' @return A named list: file name -> lines.
.kit_syntax <- function(cases) {
  do.call(c, unname(lapply(cases, function(case) case$files)))
}

#' The batch file that estimates every model in order
#' @param syntax From [.kit_syntax()].
#' @return The batch file's lines.
.kit_batch <- function(syntax) {
  models <- names(syntax)
  c("@echo off",
    "rem Estimate every model of the multilpa comparison kit, in order.",
    "rem",
    "rem Set LG to your Latent GOLD 6.1 executable. The batch flags below (/b for",
    "rem batch, /o for the output listing) are the documented command-line form;",
    "rem check them against the manual for your build before the first run.",
    "rem Step-one models must run before the Step-3 models that read their output.",
    "",
    "set LG=\"C:\\Program Files\\LatentGOLD6.1\\lg61.exe\"",
    "cd /d \"%~dp0\"",
    "",
    as.vector(rbind(sprintf("echo Estimating %s", models),
                    sprintf("%%LG%% %s /b /o %s", models, sub("[.]lgs$", ".lst", models)))),
    "",
    "echo Done. Copy every .lst and *_posteriors.txt file back to validation\\latentgold\\returned\\")
}

#' The kit's manifest: version, run order and an md5 for every file
#' @param kit The kit directory.
#' @param syntax From [.kit_syntax()].
#' @return The manifest's lines.
.kit_manifest <- function(kit, syntax) {
  files <- setdiff(sort(list.files(kit)), "MANIFEST.txt")
  c(sprintf("multilpa %s, R %s, built %s", utils::packageVersion("multilpa"),
            getRversion(), format(Sys.Date())),
    "",
    "Run order:",
    sprintf("  %d. %s", seq_along(syntax), names(syntax)),
    "",
    "md5:",
    sprintf("  %s  %s", tools::md5sum(file.path(kit, files)), files))
}

#' One row per case: what was fitted and how firmly
#' @param targets A list of targets.
#' @return A data frame.
.kit_summary <- function(targets) {
  do.call(rbind, lapply(targets, function(target) {
    data.frame(case = target$case, log_likelihood = target$log_likelihood,
               n_parameters = target$n_parameters,
               best_replicated = target$n_best_replicated,
               syntax = paste(target$lgs, collapse = " "), row.names = NULL)
  }))
}

cases <- lg_cases()
message(sprintf("Fitting %d cases with 50 starts each.", length(cases)))
targets <- Map(function(name, case) {
  message("  ", name)
  .targets(name, case, case$fit())
}, names(cases), cases)

.write_kit_data(cases, kit)
syntax <- .kit_syntax(cases)
invisible(Map(function(name, lines) .write_windows(lines, file.path(kit, name)),
              names(syntax), syntax))
invisible(Map(function(name, target) {
  saveRDS(target, file.path(targets_dir, sprintf("%s.rds", name)))
}, names(targets), targets))
.write_windows(.kit_batch(syntax), file.path(kit, "run-all.bat"))
.write_windows(.kit_manifest(kit, syntax), file.path(kit, "MANIFEST.txt"))

print(.kit_summary(targets), row.names = FALSE, digits = 10)
message(sprintf("Kit written to %s (%d files); targets to %s.", kit,
                length(list.files(kit)), targets_dir))
