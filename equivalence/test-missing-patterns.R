# Five missingness patterns, each fitted by genuine Mplus 9 (Demo) and Latent
# GOLD 6.1, against this package. The data, inputs and returned outputs are in
# missing-patterns/; see its README for the conventions the programs differ on.

missing_pattern_cases <- readRDS(file.path("missing-patterns", "kit", "cases.rds"))
missing_pattern_returned <- file.path("missing-patterns", "returned")
missing_pattern_vars <- c("y1", "y2", "y3")

fit_missing_pattern <- function(case) {
  if (identical(case$level, "single")) {
    return(lpa(case$data, missing_pattern_vars, n_profiles = case$n_profiles,
               model = "VVI", missing = "fiml", n_starts = 50, seed = 1,
               tol = 1e-10, max_iter = 5000))
  }
  multilpa(case$data, missing_pattern_vars, "g", case$n_profiles,
           case$n_group_classes, model = "VVI", missing = "fiml", n_starts = 50,
           seed = 1, tol = 1e-10, max_iter = 5000)
}

#' The number printed at the end of the first line matching `pattern`
printed_value <- function(lines, pattern) {
  line <- grep(pattern, lines, value = TRUE)[1L]
  as.numeric(utils::tail(strsplit(trimws(line), "[[:space:]]+")[[1L]], 1L))
}

#' The relabelling of this package's classes that brings its posteriors
#' closest to the reference's, and the largest difference under it
#'
#' One permutation for the whole matrix, so a row cannot be matched to a
#' different class from the rest. The same permutation then matches AvePP, so
#' the two checks agree on which class is which.
#' @return A list: `order` (this package's class for each reference class) and
#'   `gap`.
relabel <- function(mine, theirs) {
  permutations <- function(k) if (k == 1L) list(1L) else
    do.call(c, lapply(seq_len(k), \(i) lapply(permutations(k - 1L),
                                             \(p) c(i, setdiff(seq_len(k), i)[p]))))
  candidates <- permutations(ncol(mine))
  gaps <- vapply(candidates, \(p) max(abs(mine[, p, drop = FALSE] - theirs)), numeric(1))
  list(order = candidates[[which.min(gaps)]], gap = min(gaps))
}

#' Rows of the package's data a reference's ids select, after checking that
#' the ids are unique and exactly the set expected
reference_rows <- function(ids, data, expected_ids) {
  expect_false(anyDuplicated(ids) > 0L)
  expect_setequal(ids, expected_ids)
  match(ids, data$id)
}

#' Mplus's saved individual posteriors, by id; two-level joints marginalized
mplus_saved_posteriors <- function(name, case) {
  saved <- utils::read.table(file.path(missing_pattern_returned,
                                       paste0(name, "-posteriors.dat")),
                             na.strings = "*")
  if (identical(case$level, "single")) {
    return(list(id = saved[[ncol(saved)]],
                posteriors = as.matrix(saved[, 3L + seq_len(case$n_profiles)])))
  }
  # Columns: y1 y2 y3, the four cb/cw joints (11 12 21 22), cb, cw, joint, id, g.
  joint <- as.matrix(saved[, 4:7])
  list(id = saved[[ncol(saved) - 1L]],
       posteriors = cbind(joint[, 1] + joint[, 3], joint[, 2] + joint[, 4]))
}

test_that("the five patterns' data are what the kit describes", {
  observed <- lapply(missing_pattern_cases, \(case)
    rowSums(!is.na(case$data[missing_pattern_vars])))
  expect_identical(vapply(observed, \(o) sum(o == 0L), integer(1)),
                   c(m1 = 0L, m2 = 62L, m3 = 0L, m4 = 0L, m5 = 110L))
  m5 <- missing_pattern_cases$m5$data
  expect_identical(sort(unique(m5$g[rowSums(!is.na(m5[missing_pattern_vars])) > 0])),
                   setdiff(1:60, c(7L, 31L)))
})

test_that("likelihood, sample size, BIC, entropy and AvePP reproduce Mplus", {
  invisible(lapply(names(missing_pattern_cases), function(name) {
    case <- missing_pattern_cases[[name]]
    fit <- fit_missing_pattern(case)
    path <- file.path(missing_pattern_returned, paste0(name, ".out"))
    lines <- read_mplus_lines(path)
    expect_true(any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)))
    expect_true(any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)))
    # Mplus prints three decimals; two searches stop at their own criteria.
    expect_lt(abs(fit$log_likelihood - printed_value(lines, "H0 Value")),
              5e-4 + 1e-5, label = name)
    # Mplus leaves out rows with nothing observed, as this package's sample
    # size does.
    expect_identical(fit$n_informative,
                     as.integer(printed_value(lines, "^Number of observations")))
    expect_lt(abs(fit$bic_individual - printed_value(lines, "Bayesian \\(BIC\\)")),
              5e-4 + 1e-4, label = name)
    theirs <- if (identical(case$level, "single"))
      mplus_classification_single(path, case$n_profiles) else
        mplus_classification_twolevel(path, case$n_profiles,
                                      case$n_group_classes)$individuals
    # Mplus saves exactly the rows it analysed: those with an observed value.
    saved <- mplus_saved_posteriors(name, case)
    observed <- rowSums(!is.na(case$data[missing_pattern_vars])) > 0L
    rows <- reference_rows(saved$id, case$data, case$data$id[observed])
    matched <- relabel(fit$subject_posteriors[rows, , drop = FALSE], saved$posteriors)
    expect_lt(matched$gap, 5e-4, label = paste(name, "posteriors"))
    # Entropy and AvePP, with classes matched by the posteriors above.
    discrepancy <- classification_discrepancy(
      theirs, latents_classification(fit, "individuals"), matched$order)
    expect_lt(discrepancy[["entropy"]], mplus_printed_tolerance,
              label = paste(name, "entropy"))
    expect_lt(discrepancy[["avepp"]], mplus_printed_tolerance,
              label = paste(name, "AvePP"))
  }))
})

test_that("likelihood and every row's posteriors reproduce Latent GOLD", {
  invisible(lapply(names(missing_pattern_cases), function(name) {
    case <- missing_pattern_cases[[name]]
    fit <- fit_missing_pattern(case)
    lines <- iconv(readLines(file.path(missing_pattern_returned, paste0(name, ".lst")),
                             warn = FALSE), "latin1", "UTF-8")
    # Latent GOLD prints four decimals.
    expect_lt(abs(fit$log_likelihood - printed_value(lines, "^Log-likelihood \\(LL\\)")),
              5e-5 + 1e-5, label = name)
    saved <- utils::read.delim(file.path(missing_pattern_returned,
                                         paste0(name, "_posteriors.txt")),
                               check.names = FALSE)
    theirs <- as.matrix(saved[, paste0("Cluster#", seq_len(case$n_profiles))])
    # Latent GOLD keeps every row, including those with nothing observed, and
    # gives them, as this package does, the prior their group implies.
    rows <- reference_rows(saved$id, case$data, case$data$id)
    expect_lt(relabel(fit$subject_posteriors[rows, , drop = FALSE], theirs)$gap,
              5e-4, label = paste(name, "posteriors"))
  }))
})

test_that("class counts reproduce Mplus's counts based on posterior probabilities", {
  # Mplus counts the classified rows only; so does drop_empty = TRUE.
  invisible(lapply(names(missing_pattern_cases), function(name) {
    case <- missing_pattern_cases[[name]]
    fit <- fit_missing_pattern(case)
    lines <- read_mplus_lines(file.path(missing_pattern_returned, paste0(name, ".out")))
    counts <- if (identical(case$level, "single"))
      mplus_posterior_counts(lines, case$n_profiles) else
        mplus_variable_counts(lines, "CW", case$n_profiles)
    saved <- mplus_saved_posteriors(name, case)
    rows <- match(saved$id, case$data$id)
    matched <- relabel(fit$subject_posteriors[rows, , drop = FALSE], saved$posteriors)
    # Mplus prints the counts to five decimals.
    expect_lt(max(abs(fit$effective_profile_counts[matched$order] - counts)), 0.05,
              label = name)
  }))
})

test_that("the programs' sample sizes differ only by the rows with nothing observed", {
  invisible(lapply(c("m2", "m5"), function(name) {
    case <- missing_pattern_cases[[name]]
    empty <- sum(rowSums(!is.na(case$data[missing_pattern_vars])) == 0L)
    mplus <- read_mplus_lines(file.path(missing_pattern_returned, paste0(name, ".out")))
    latent_gold <- iconv(readLines(file.path(missing_pattern_returned, paste0(name, ".lst")),
                                   warn = FALSE), "latin1", "UTF-8")
    expect_identical(as.integer(printed_value(mplus, "^Number of observations")),
                     nrow(case$data) - empty)
    expect_identical(as.integer(printed_value(latent_gold, "^Number of cases")),
                     nrow(case$data))
  }))
})

test_that("m5's group level reproduces Mplus: clusters, posteriors and statistics", {
  # Groups with nothing observed are left out, as Mplus leaves them out of its
  # clusters. Mplus then averages the group-level entropy and AvePP over
  # observed rows, each carrying its group's posterior; this package averages
  # over groups, each counted once. The Mplus figures are therefore checked by
  # applying Mplus's averaging to this package's classified group posteriors.
  case <- missing_pattern_cases$m5
  fit <- fit_missing_pattern(case)
  lines <- read_mplus_lines(file.path(missing_pattern_returned, "m5.out"))
  expect_identical(fit$n_informative_groups,
                   as.integer(printed_value(lines, "Number of clusters")))
  expect_identical(fit$n_groups, 60L)
  printed <- mplus_classification_twolevel(file.path(missing_pattern_returned, "m5.out"),
                                           2L, 2L)$groups
  per_observed_row <- fit$group_posteriors[fit$group_index[fit$informative], , drop = FALSE]
  expect_lt(abs(.multilpa_relative_entropy(per_observed_row) - printed$entropy),
            mplus_printed_tolerance)
  # The package's own group entropy is over the 58 classified groups.
  groups <- subset(get_results(fit, "entropy"), level == "groups")
  expect_identical(groups$n_units, 58L)
  observed_rows <- rowSums(!is.na(case$data[missing_pattern_vars])) > 0L
  # Mplus: the cb margin of the saved joints, one row per analysed group.
  saved <- utils::read.table(file.path(missing_pattern_returned, "m5-posteriors.dat"),
                             na.strings = "*")
  joint <- as.matrix(saved[, 4:7])
  per_row <- cbind(joint[, 1] + joint[, 2], joint[, 3] + joint[, 4])
  group_of_row <- saved[[ncol(saved)]]
  first <- !duplicated(group_of_row)
  expect_setequal(group_of_row[first], unique(case$data$g[observed_rows]))
  mine <- fit$group_posteriors[match(group_of_row[first], fit$group_values), , drop = FALSE]
  matched <- relabel(mine, per_row[first, , drop = FALSE])
  expect_lt(matched$gap, 5e-4)
  # Mplus counts the group classes over observed rows; the package counts
  # groups. Weighting the package's group posteriors by each group's observed
  # rows gives Mplus's figures.
  row_counts <- colSums(per_observed_row)
  expect_lt(max(abs(row_counts[matched$order] - mplus_variable_counts(lines, "CB", 2L))), 0.05)
  expect_equal(sum(fit$effective_group_counts), 58, tolerance = 1e-8)
  # AvePP under Mplus's averaging, classes matched by the group posteriors.
  avepp <- diag(.multilpa_average_posterior_matrix(per_observed_row))
  expect_lt(max(abs(avepp[matched$order] - printed$avepp)), mplus_printed_tolerance)
  # Latent GOLD: every group, including the two with nothing observed.
  lg <- utils::read.delim(file.path(missing_pattern_returned, "m5_posteriors.txt"),
                          check.names = FALSE)
  lg_groups <- lg[!duplicated(lg[[1L]]), , drop = FALSE]
  expect_setequal(lg_groups[[1L]], unique(case$data$g))
  mine_all <- fit$group_posteriors[match(lg_groups[[1L]], fit$group_values), , drop = FALSE]
  expect_lt(relabel(mine_all, as.matrix(lg_groups[, c("GClass#1", "GClass#2")]))$gap, 5e-4)
})
