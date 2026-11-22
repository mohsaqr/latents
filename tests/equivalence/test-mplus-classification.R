# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# Classification summaries against Mplus 9 (Muthen and Muthen), from two runs
# kept with the package's verification material (fixtures/mplus-avepp.json,
# written by fixtures/make-mplus-avepp.R from the .out files; the data beside
# it in data/): a three-class latent profile model with equal variances and no
# covariances (Mplus's default, mclust's EEI) on 500 poorly separated rows, and
# a three-class latent class model of six binary items on 1000 rows. Mplus ran
# 200 random starts (40 final stages). Compared: the H0 log likelihood, the
# relative entropy, the estimated class counts and the AvePP matrix (average
# latent class probabilities for most likely class membership).
#
# Mplus prints three decimals for the likelihood, the entropy and the AvePP:
# half a unit of the last digit, 5e-4 absolute. The class counts are printed to
# five decimals, but Mplus stops at a likelihood change of 1e-8 and a
# derivative of 1e-6, so its counts sit within about 1e-3 of the maximum:
# 5e-3 absolute. Mplus's class labels are matched to latents' by the order of
# the estimated counts that is closest to Mplus's.

mplus <- reference_fixture("fixtures/mplus-avepp.json")

mplus_data <- function(case, names) {
  data <- utils::read.table(reference_data(sprintf("data/mplus-%s.dat", case)), col.names = names)
  data[data == -999] <- NA
  data
}

all_orders <- function(k) {
  if (k == 1L) return(list(1L))
  do.call(c, lapply(seq_len(k), \(first) lapply(all_orders(k - 1L), \(rest) c(first, setdiff(seq_len(k), first)[rest]))))
}

check_with_mplus <- function(fit, reference, label, covers_fit) {
  expect_reference(as.numeric(logLik(fit)), reference$log_likelihood, reference = mplus, authority = NULL,
                   covers = covers_fit, quantity = sprintf("log likelihood, %s", label), tolerance = 5e-4,
                   measure = "absolute")
  classification <- subset(get_classification(fit), level == "individuals")
  # order[m] is latents' label of Mplus's class m.
  orders <- all_orders(3L)
  gaps <- vapply(orders, \(order) sum((classification$estimated_n[order] - reference$class_counts$estimated_n)^2),
                 numeric(1))
  order <- orders[[which.min(gaps)]]
  avepp <- subset(get_avepp(fit), level == "individuals")
  mplus_label <- match(seq_len(3L), order)
  relabelled <- data.frame(assigned_class = mplus_label[avepp$assigned_class], class = mplus_label[avepp$class],
                           average_posterior = avepp$average_posterior)
  expect_reference(relabelled, reference$avepp, reference = mplus, covers = "get_avepp()",
                   quantity = sprintf("AvePP matrix, %s", label), tolerance = 5e-4, measure = "absolute",
                   key = c("assigned_class", "class"))
  expect_reference(subset(get_entropy(fit), level == "individuals")$relative_entropy, reference$relative_entropy,
                   reference = mplus, covers = "get_entropy()", quantity = sprintf("relative entropy, %s", label),
                   tolerance = 5e-4, measure = "absolute")
  diagnosed <- diagnostics(fit, plots = FALSE)
  average <- get_results(diagnosed, "average_posteriors")
  expect_reference(data.frame(assigned_class = mplus_label[average$assigned_class], class = mplus_label[average$class],
                              average_posterior = average$average_posterior),
                   reference$avepp, reference = mplus, covers = c("diagnostics()", "get_results.multilpa_diagnostics()"),
                   quantity = sprintf("AvePP matrix, diagnostics(), %s", label), tolerance = 5e-4,
                   measure = "absolute", key = c("assigned_class", "class"))
  # latents' estimated class counts in Mplus's class order, for the caller to compare.
  data.frame(class = 1:3, estimated_n = classification$estimated_n[order])
}

test_that("AvePP, entropy and class counts of a three-profile EEI model are Mplus's", {
  data <- mplus_data("eei_k3poor", c("id", "y1", "y2", "y3", "y4"))
  fit <- lpa(data, c("y1", "y2", "y3", "y4"), n_profiles = 3, model = "EEI", n_starts = 50, seed = 1,
             tol = 1e-12, max_iter = 100000L)
  counts <- check_with_mplus(fit, mplus$eei_k3poor, "three profiles, EEI, 500 rows", c("lpa()", "logLik.multilpa()"))
  expect_reference(counts, mplus$eei_k3poor$class_counts[c("class", "estimated_n")], reference = mplus,
                   covers = "get_classification()", quantity = "estimated class counts, three profiles, EEI, 500 rows",
                   tolerance = 5e-3, measure = "absolute", key = "class")
})

test_that("AvePP, entropy and class counts of a three-class model of binary items are Mplus's", {
  data <- mplus_data("lca_k3", c("id", sprintf("u%d", 1:6)))
  fit <- lca(data, sprintf("u%d", 1:6), n_classes = 3, n_starts = 50, seed = 1, tol = 1e-12, max_iter = 100000L)
  # The class counts of this model sit on a ridge of the likelihood: latents
  # (tol 1e-12) gives 373.306 / 128.179 / 498.515, poLCA (tol 1e-12, 30 starts)
  # 373.318 / 128.170 / 498.513 at a likelihood 1.4e-7 lower, Mplus 373.329 /
  # 128.160 / 498.511 (stopping at a likelihood change of 1e-8), and latents
  # stopped at a relative change of 1e-6 moves class 1 by 0.08 at a likelihood
  # only 3.7e-4 lower.
  # Measured after the first run, where the 5e-3 bound proved tighter than the
  # ridge allows; recorded as a known difference, not loosened.
  counts <- check_with_mplus(fit, mplus$lca_k3, "three classes, six binary items, 1000 rows",
                             c("lca()", "logLik.multilpa()"))
  expect_known_difference(counts, mplus$lca_k3$class_counts[c("class", "estimated_n")], reference = mplus,
                          covers = "get_classification()",
                          quantity = "estimated class counts, three classes, six binary items, 1000 rows",
                          tolerance = 5e-3, measure = "absolute", key = "class",
                          reason = paste(
                            "The counts lie on a ridge of the likelihood flat to 1e-7: latents, poLCA and Mplus stop",
                            "within 0.02 of each other (373.306, 373.318, 373.329 for class 1) with the same",
                            "likelihood to Mplus's three printed decimals; Mplus stops at a likelihood change of 1e-8."))
})
