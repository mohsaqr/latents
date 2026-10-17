# Reading Mplus's classification-quality output: the relative entropy and the
# average latent class probabilities for most likely class membership (AvePP),
# which Mplus prints to three decimals. The comparison is with this package's
# own get_results() tables, not with a formula written in the test.

#' The lines of an Mplus output, with any HTML markup of a published copy removed
#' @param path An Mplus `.out` file or a published `.html` copy of one.
#' @return A character vector of lines.
read_mplus_lines <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (grepl("[.]html$", path)) lines <- gsub("<[^>]+>", "", lines)
  gsub("&nbsp;", " ", lines, fixed = TRUE)
}

#' The first `Entropy` value printed after line `from`
#' @param lines Output lines.
#' @param from Line index to search after.
#' @return A single number.
mplus_entropy_after <- function(lines, from) {
  at <- grep("^\\s*Entropy\\s+[0-9.]+\\s*$", lines)
  at <- at[at > from][1L]
  stopifnot("no Entropy line after the anchor" = !is.na(at))
  as.numeric(sub(".*Entropy\\s+", "", lines[at]))
}

#' The diagonal of the first AvePP matrix printed after line `from`
#'
#' Mplus prints a column header row (`1 2 ...`) and then one row per most
#' likely class, each starting with the class label. Row `k`, column `k` is the
#' AvePP of class `k`, in Mplus's own class labels.
#' @param lines Output lines.
#' @param from Line index to search after.
#' @param k Number of classes.
#' @return A numeric vector of length `k`.
mplus_avepp_after <- function(lines, from, k) {
  header <- grep("Average Latent Class Probabilities for Most Likely", lines)
  start <- header[header > from][1L]
  stopifnot("no AvePP table after the anchor" = !is.na(start))
  numeric_rows <- grep("^\\s+[0-9]+(\\s+-?[0-9.]+)+\\s*$", lines)
  rows <- numeric_rows[numeric_rows > start][seq_len(k + 1L)]
  values <- do.call(rbind, lapply(lines[rows[-1L]], \(line)
    as.numeric(strsplit(trimws(line), "\\s+")[[1L]])[-1L]))
  stopifnot("the AvePP table is not k by k" = identical(dim(values), c(k, k)))
  diag(values)
}

#' Class counts based on estimated posterior probabilities, in Mplus's labels
#'
#' These fix the correspondence between Mplus's class labels and this
#' package's, which may number the same classes differently.
#' @param lines Output lines.
#' @param k Number of classes.
#' @return A numeric vector of length `k`.
mplus_posterior_counts <- function(lines, k) {
  at <- grep("BASED ON ESTIMATED POSTERIOR PROBABILITIES", lines)[1L]
  rows <- grep("^\\s+[0-9]+\\s+[0-9.]+\\s+[0-9.]+\\s*$", lines)
  rows <- rows[rows > at][seq_len(k)]
  vapply(lines[rows], \(line) as.numeric(strsplit(trimws(line), "\\s+")[[1L]][2L]),
         numeric(1), USE.NAMES = FALSE)
}

#' Classification quality of a single-level Mplus output
#' @param path The output.
#' @param k Number of classes.
#' @return A list with `entropy`, `avepp` and `counts`.
mplus_classification_single <- function(path, k) {
  lines <- read_mplus_lines(path)
  quality <- grep("CLASSIFICATION QUALITY", lines)[1L]
  list(entropy = mplus_entropy_after(lines, quality),
       avepp = mplus_avepp_after(lines, quality, k),
       counts = mplus_posterior_counts(lines, k))
}

#' Classification quality of each level of a two-level Mplus output
#'
#' Mplus names the within (profile) variable `CW` and the between (group
#' class) variable `CB`, each with its own entropy and AvePP table.
#' @param path The output.
#' @param k_profiles,k_groups Numbers of classes at each level.
#' @return A list with `individuals` and `groups`, each `entropy` and `avepp`.
mplus_classification_twolevel <- function(path, k_profiles, k_groups) {
  lines <- read_mplus_lines(path)
  within <- grep("Classification Quality for CW\\s*$", lines)[1L]
  between <- grep("Classification Quality for CB\\s*$", lines)[1L]
  stopifnot("no CW or CB classification block" = !is.na(within) && !is.na(between))
  list(individuals = list(entropy = mplus_entropy_after(lines, within),
                          avepp = mplus_avepp_after(lines, within, k_profiles)),
       groups = list(entropy = mplus_entropy_after(lines, between),
                     avepp = mplus_avepp_after(lines, between, k_groups)))
}

#' This package's reported entropy and AvePP at one level
#' @param fit A fitted model.
#' @param which `"individuals"` or `"groups"`.
#' @return A list with `entropy` and `avepp`, AvePP in the fit's class order.
latents_classification <- function(fit, which) {
  entropy <- subset(get_results(fit, "entropy"), level == which)
  classes <- subset(get_results(fit, "classification"), level == which)
  list(entropy = entropy$relative_entropy, avepp = classes$average_posterior)
}

#' Largest discrepancy between Mplus's printed values and this package's
#' @param theirs,mine Lists with `entropy` and `avepp`.
#' @param order This package's class index for each Mplus class.
#' @return A named numeric vector: `entropy` and `avepp`.
classification_discrepancy <- function(theirs, mine, order) {
  c(entropy = abs(theirs$entropy - mine$entropy),
    avepp = max(abs(theirs$avepp - mine$avepp[order])))
}

# Mplus prints three decimals, so a value it prints agrees with any value
# within half a unit of the last place.
mplus_printed_tolerance <- 5e-4
