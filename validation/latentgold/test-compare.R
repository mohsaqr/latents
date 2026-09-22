# Tests for the Latent GOLD comparison, run before any real output exists.
#
# Latent GOLD output is imitated from this package's own targets: posteriors
# rounded to four decimals and written with the class labels swapped, and a
# listing holding the log likelihood and the parameter count. The comparison
# must call that agreement, and must call each deliberate break a disagreement
# or say it could not read the file. A comparison that only ever says "agree"
# would pass the first test and fail every other one.
#
# From the project root, after make-kit.R:
#   Rscript validation/latentgold/test-compare.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("validation", "latentgold", "compare-functions.R"))
library(testthat)

targets_dir <- file.path("validation", "latentgold", "targets")
if (!file.exists(file.path(targets_dir, "c01_continuous_varying.rds"))) {
  stop("Run validation/latentgold/make-kit.R first; these tests read its targets.")
}
target <- readRDS(file.path(targets_dir, "c01_continuous_varying.rds"))
sequences <- readRDS(file.path(targets_dir, "c11_lta_mixture.rds"))

#' Imitate what Latent GOLD would write for one case
#'
#' @param target One case's targets.
#' @param directory Where to write.
#' @param individual,group Latent variable names for the posterior columns.
#' @param digits Decimals to round the posteriors to.
#' @param separator Field separator of the posterior file.
#' @param likelihood,parameters Values for the listing.
#' @return `directory`, invisibly.
imitate <- function(target, directory, individual = "Cluster", group = "GClass",
                    digits = 4L, separator = "\t",
                    likelihood = round(target$log_likelihood, 4L),
                    parameters = target$n_parameters,
                    n_cases = target$n_observations,
                    n_groups = if (!is.null(target$group_posteriors)) target$n_groups) {
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  stem <- target$case
  # Latent GOLD writes each criterion once per sample size it uses: the
  # unsuffixed family on its cases, and the ",Ngroups" family only when the
  # model has a group level.
  criterion <- function(label, column, suffix) {
    value <- target$criteria[[column]]
    if (is.null(value)) NULL else sprintf("%s (based on LL%s)\t%.4f\t\t\t", label, suffix, value)
  }
  convention <- if (isTRUE(all.equal(as.numeric(n_cases), as.numeric(target$n_observations))))
    "individual" else "groups"
  criteria_lines <- c(
    criterion("AIC", "aic", ""),
    criterion("BIC", paste0("bic_", convention), ""),
    criterion("CAIC", paste0("caic_", convention), ""),
    criterion("SABIC", paste0("sabic_", convention), ""),
    if (!is.null(n_groups)) c(criterion("BIC", "bic_groups", ",Ngroups"),
                              criterion("CAIC", "caic_groups", ",Ngroups"),
                              criterion("SABIC", "sabic_groups", ",Ngroups")))
  # The labels, the tab separator and the Latin-1 "R\u00b2" are copied from a real
  # Latent GOLD 6.1 listing produced under Wine on 2026-09-22. "Npar" also
  # appears alone as a column heading of the summary table, on a line with no
  # value: the parser must not read that one.
  listing <- c("c01.dat",
               "\t\tLL\tBIC(LL)\tAIC(LL)\tNpar\tClass.Err.\tEntropy R\xb2\t",
               sprintf("model\t2-GClass  2-Cluster\t%.4f\t0\t0\t%d\t0.11\t0.62\t",
                       likelihood, parameters),
               "Log-likelihood Statistics\t\t\t\t",
               sprintf("Log-likelihood (LL)\t%.4f\t\t\t", likelihood),
               "Log-prior\t0.0000\t\t\t",
               sprintf("Number of parameters (Npar)\t%d\t\t\t", parameters),
               sprintf("Number of cases\t%d\t\t\t", n_cases),
               if (!is.null(n_groups)) sprintf("Number of groups\t%d\t\t\t", n_groups),
               criteria_lines)
  connection <- file(file.path(directory, sprintf("%s.lst", stem)), open = "wb")
  writeLines(listing, connection, useBytes = TRUE)
  close(connection)
  classes <- as.matrix(target$posteriors[grep("^class_", names(target$posteriors))])
  # Latent GOLD's labels are its own: reverse them so alignment is exercised.
  posterior <- formatC(classes[, rev(seq_len(ncol(classes))), drop = FALSE],
                       digits = digits, format = "f")
  frame <- data.frame(id = target$posteriors$id, posterior, check.names = FALSE)
  names(frame) <- c("id", sprintf("%s#%d", individual, seq_len(ncol(classes))))
  if (!is.null(target$group_posteriors)) {
    group_of_row <- target$data$g[match(target$posteriors$id, target$data$id)]
    groups <- as.matrix(target$group_posteriors[grep("^class_", names(target$group_posteriors))])
    rows <- match(group_of_row, target$group_posteriors$g)
    labelled <- formatC(groups[rows, rev(seq_len(ncol(groups))), drop = FALSE],
                        digits = digits, format = "f")
    colnames(labelled) <- sprintf("%s#%d", group, seq_len(ncol(groups)))
    frame <- cbind(frame, labelled)
  }
  utils::write.table(frame, file.path(directory, sprintf("%s_posteriors.txt", stem)),
                     sep = separator, quote = FALSE, row.names = FALSE)
  invisible(directory)
}

statuses_of <- function(result, pattern) {
  result$status[grepl(pattern, result$quantity)]
}

# --------------------------------------------------------------------------
test_that("small helpers do what the comparison relies on", {
  expect_equal(.lg_half_unit(c("0.1234", "-12.5")), 5e-5)
  # A value printed without its trailing zeros does not coarsen the file.
  expect_equal(.lg_half_unit(c("0.1234", "1", "0")), 5e-5)
  expect_equal(.lg_half_unit("-2790.0958"), 5e-5)
  expect_equal(.lg_half_unit("15"), 0.5)
  orders <- .lg_permutations(3L)
  expect_equal(nrow(unique(orders)), 6L)
  expect_true(all(apply(orders, 1L, function(row) setequal(row, 1:3))))
  swapped <- .lg_align(diag(2), diag(2)[, 2:1])
  expect_equal(swapped$order, 2:1)
  expect_equal(swapped$difference, 0)
})

test_that("imitated Latent GOLD output agrees on every compared quantity", {
  returned <- imitate(target, tempfile("lg-agree-"))
  result <- lg_compare_case(target, returned)
  compared <- result[result$status %in% c("agree", "disagree"), , drop = FALSE]
  # likelihood, parameter count, 4 + 3 information criteria, two sample sizes,
  # two posterior levels, 2 x 3 means, 2 x 3 variances, two group proportions
  expect_equal(nrow(compared), 27L)
  expect_true(all(compared$status == "agree"),
              info = paste(compared$quantity[compared$status != "agree"], collapse = "; "))
  expect_true(all(result$status[!result$status %in% c("agree", "disagree")] ==
                    "awaiting parser"))
})

test_that("the moment tolerance is derived from the rounding, and holds at it", {
  # Rounding to two decimals is a far coarser file. The derived bound widens
  # with it, and the recomputed means and variances still fall inside.
  returned <- imitate(target, tempfile("lg-coarse-"), digits = 2L)
  result <- lg_compare_case(target, returned)
  moments <- result[grepl("^(mean|variance):", result$quantity), , drop = FALSE]
  fine <- lg_compare_case(target, imitate(target, tempfile("lg-fine-")))
  fine_moments <- fine[grepl("^(mean|variance):", fine$quantity), , drop = FALSE]
  expect_true(all(moments$status == "agree"))
  expect_true(all(moments$tolerance > fine_moments$tolerance))
})

test_that("a comma-separated posterior file is read the same way", {
  returned <- imitate(target, tempfile("lg-comma-"), separator = ",")
  result <- lg_compare_case(target, returned)
  expect_true(all(statuses_of(result, "^posteriors") == "agree"))
})

test_that("a wrong likelihood or parameter count is a disagreement", {
  returned <- imitate(target, tempfile("lg-wrong-"),
                      likelihood = round(target$log_likelihood, 4L) - 0.01,
                      parameters = target$n_parameters + 1L)
  result <- lg_compare_case(target, returned)
  expect_equal(statuses_of(result, "^log likelihood$"), "disagree")
  expect_equal(statuses_of(result, "^parameters$"), "disagree")
})

test_that("a posterior moved by 0.01 is a disagreement", {
  returned <- imitate(target, tempfile("lg-moved-"))
  path <- file.path(returned, sprintf("%s_posteriors.txt", target$case))
  # `#` starts a comment by default, and every posterior column name holds one.
  frame <- utils::read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                             comment.char = "")
  frame[["Cluster#1"]][[1L]] <- frame[["Cluster#1"]][[1L]] + 0.01
  frame[["Cluster#2"]][[1L]] <- frame[["Cluster#2"]][[1L]] - 0.01
  utils::write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE)
  result <- lg_compare_case(target, returned)
  expect_equal(statuses_of(result, "^posteriors: individuals"), "disagree")
  expect_equal(statuses_of(result, "^posteriors: groups"), "agree")
})

test_that("absent output is missing, and unreadable output says what it found", {
  empty <- tempfile("lg-empty-")
  dir.create(empty)
  absent <- lg_compare_case(target, empty)
  expect_true(all(absent$status %in% c("missing", "awaiting parser")))
  expect_true(any(absent$status == "missing"))

  renamed <- imitate(target, tempfile("lg-renamed-"), individual = "Clu")
  result <- lg_compare_case(target, renamed)
  unread <- result[grepl("^posteriors: individuals", result$quantity), , drop = FALSE]
  expect_equal(unread$status, "unparsed")
  expect_match(unread$note, "Clu#1")
})

test_that("the listing is read as Latin-1, and Npar's column heading is not a value", {
  # Reading the listing as UTF-8 leaves invalid strings that match nothing, and
  # the summary table's "Npar" heading carries no value of its own.
  returned <- imitate(target, tempfile("lg-latin1-"))
  result <- lg_compare_case(target, returned)
  expect_equal(statuses_of(result, "^log likelihood$"), "agree")
  expect_equal(statuses_of(result, "^parameters$"), "agree")
  lines <- .lg_listing_lines(file.path(returned, sprintf("%s.lst", target$case)))
  expect_true(any(grepl("Entropy R", lines)))
  expect_false(any(is.na(nchar(lines, allowNA = TRUE))))
})

test_that("a listing with no summary lines is unparsed, not agreed", {
  returned <- imitate(target, tempfile("lg-blank-"))
  writeLines("Latent GOLD crashed", file.path(returned, sprintf("%s.lst", target$case)))
  result <- lg_compare_case(target, returned)
  expect_equal(statuses_of(result, "^log likelihood$"), "unparsed")
  expect_equal(statuses_of(result, "^parameters$"), "unparsed")
})

test_that("a latent transition case reads State and GClass posteriors", {
  returned <- imitate(sequences, tempfile("lg-lta-"), individual = "State")
  result <- lg_compare_case(sequences, returned)
  expect_true(all(statuses_of(result, "^posteriors") == "agree"))
  expect_equal(length(statuses_of(result, "^posteriors")), 2L)
  expect_true(all(statuses_of(result, "^(mean|variance):") == "agree"))
})

test_that("a criterion is matched to the column computed on the same N", {
  # Latent GOLD's "case" is an observation in a two-level model but a whole
  # sequence in a caseid model, which is this package's group. Matching by
  # label instead of by sample size compared a sequence model's BIC against
  # bic_individual and called it a disagreement.
  two_level <- lg_compare_case(target, imitate(target, tempfile("lg-n-")))
  expect_equal(statuses_of(two_level, "^(bic|caic|sabic)_(individual|groups)$"),
               rep("agree", 6L))
  expect_equal(statuses_of(two_level, "^n_(observations|groups)$"), c("agree", "agree"))

  # A sequence model: Latent GOLD counts sequences and prints no group family.
  returned <- imitate(sequences, tempfile("lg-seq-"), individual = "State",
                      n_cases = sequences$n_groups, n_groups = NULL)
  result <- lg_compare_case(sequences, returned)
  expect_equal(statuses_of(result, "^(bic|caic|sabic)_groups$"), rep("agree", 3L))
  expect_false(any(grepl("_individual$", result$quantity)))
  expect_equal(statuses_of(result, "^n_groups$"), "agree")
})

test_that("output produced from a different kit is reported as stale", {
  kit <- tempfile("lg-kit-"); returned <- tempfile("lg-ret-")
  dir.create(kit); dir.create(returned)
  writeLines("id\ty1", file.path(kit, "c01.dat"))
  writeLines("infile 'c01.dat'", file.path(kit, "c01.lgs"))
  writeLines(sprintf("%s  %s", unname(lg_fingerprint(kit)), names(lg_fingerprint(kit))),
             file.path(returned, "INPUTS.md5"))
  expect_true(lg_check_fingerprint(kit, returned))
  # The data change; the output in returned/ no longer describes this kit.
  writeLines("id\ty1\ty2", file.path(kit, "c01.dat"))
  expect_warning(result <- lg_check_fingerprint(kit, returned),
                 class = "multilpa_stale_latentgold_output")
  expect_false(result)
  # A kit gaining a case is a change too.
  expect_warning(lg_check_fingerprint(kit, tempfile("lg-none-")), NA)
})

test_that("every generated syntax line is well formed", {
  # Latent GOLD rejected `GClass <- 1 + ;` in seven cases on 2026-09-22:
  # paste0(" + ", character(0), collapse = "") is " + ", not "". Only the one
  # case that had covariates escaped, so a check of one file would have missed
  # it; every file is checked here.
  source(file.path("validation", "latentgold", "cases.R"), local = TRUE)
  syntax <- do.call(c, unname(lapply(lg_cases(), function(case) case$files)))
  expect_gt(length(syntax), 10L)
  for (name in names(syntax)) {
    lines <- syntax[[name]]
    dangling <- grep("[+|,]\\s*;\\s*$|<-\\s*;", lines, value = TRUE)
    expect_equal(dangling, character(), info = sprintf("%s: %s", name,
                                                       paste(dangling, collapse = " / ")))
    expect_equal(sum(lines == "model"), sum(lines == "end model"),
                 info = name)
    expect_true(any(grepl("^infile ", lines)), info = name)
  }
})

cat("\nLatent GOLD comparison self-tests complete.\n")
