# enumerate_lpa() on the Old Faithful eruptions (datasets::faithful) with one to
# four profiles under mclust's EEI and VVI structures, and the accessors that
# read the grid. Reference: mclust::mclustBIC() on the same grid, both programs
# run to a relative change of 1e-12; the criteria, the differences between
# consecutive candidates, the best candidate per criterion and the knee of the
# BIC line are then computed from mclust's numbers by their definitions.

enumeration_fit <- enumerate_lpa(datasets::faithful, c("eruptions", "waiting"), n_profiles = 1:4,
                                 model = c("EEI", "VVI"), n_starts = 20, seed = 1, tol = 1e-12,
                                 max_iter = 100000L)
enumeration_n <- nrow(datasets::faithful)

# mclust's grid as one row per candidate, named as latents names its candidates.
mclust_grid <- function() {
  bic <- mclust::mclustBIC(datasets::faithful, G = 1:4, modelNames = c("EEI", "VVI"),
                           control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(1e6, 1e6)),
                           verbose = FALSE)
  cells <- expand.grid(n_profiles = 1:4, model = c("EEI", "VVI"), stringsAsFactors = FALSE)
  cells$candidate <- sprintf("profiles%d_%s", cells$n_profiles, cells$model)
  cells$n_parameters <- mapply(\(model, k) mclust::nMclustParams(model, 2L, k), cells$model, cells$n_profiles)
  # mclust's BIC is 2 log L - p log n: its sign is the opposite of latents'.
  cells$bic_individual <- -mapply(\(model, k) bic[as.character(k), model], cells$model, cells$n_profiles)
  cells$log_likelihood <- (-cells$bic_individual + cells$n_parameters * log(enumeration_n)) / 2
  cells
}

# The candidate where mclust's hierarchical initialisation stops at a lower maximum (see the known difference).
vvi_three <- "profiles3_VVI"

test_that("the grid's likelihoods, parameter counts and BIC equal mclust::mclustBIC()", {
  skip_without_reference("mclust")
  grid <- mclust_grid()
  agreed <- subset(grid, candidate != vvi_three, select = c(candidate, n_parameters, log_likelihood, bic_individual))
  # Two EM runs stopped at a relative change of 1e-12: the maxima agree to about 1e-10 relative.
  expect_reference(subset(as.data.frame(enumeration_fit), candidate != vvi_three), agreed, reference = "mclust",
                   covers = c("enumerate_lpa()", "as.data.frame.multilpa_enumeration()"),
                   quantity = "log likelihood, parameters and BIC per candidate", tolerance = 1e-8,
                   measure = "relative", key = "candidate")
  expect_reference(subset(get_candidates(enumeration_fit), candidate != vvi_three), agreed, reference = "mclust",
                   covers = c("enumerate_lpa()", "get_candidates()"),
                   quantity = "log likelihood, parameters and BIC per candidate, get_candidates()",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  expect_reference(subset(get_results(enumeration_fit), candidate != vvi_three), agreed, reference = "mclust",
                   covers = c("get_results.multilpa_enumeration()", "get_results()"),
                   quantity = "log likelihood, parameters and BIC per candidate, get_results()",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  expect_reference(subset(as.data.frame(summary(enumeration_fit)), candidate != vvi_three), agreed, reference = "mclust",
                   covers = c("summary.multilpa_enumeration()", "as.data.frame.summary_multilpa_enumeration()"),
                   quantity = "log likelihood, parameters and BIC per candidate, summary()",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  expect_reference(subset(get_results(summary(enumeration_fit)), candidate != vvi_three), agreed, reference = "mclust",
                   covers = "get_results.summary_multilpa_enumeration()",
                   quantity = "log likelihood, parameters and BIC per candidate, get_results(summary())",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  expect_known_difference(subset(as.data.frame(enumeration_fit), candidate == vvi_three, select = c(candidate, log_likelihood)),
                          subset(grid, candidate == vvi_three, select = c(candidate, log_likelihood)),
                          reference = "mclust", covers = "enumerate_lpa()",
                          quantity = "log likelihood, three VVI profiles", tolerance = 1e-8, measure = "relative",
                          key = "candidate",
                          reason = "mclustBIC() starts EM from one hierarchical-clustering partition and stops at a lower local maximum for three VVI profiles; latents' random starts reach a higher one, which mclust's own EM keeps when started from latents' posteriors (next test).")
})

test_that("three VVI profiles: mclust's EM started from latents' posteriors stays at latents' maximum", {
  skip_without_reference("mclust")
  candidate <- get_candidate(enumeration_fit, n_profiles = 3, model = "VVI")
  z <- matrix(get_results(candidate, "posteriors")$posterior, ncol = 3)
  continued <- mclust::meVVI(data = as.matrix(datasets::faithful), z = z,
                             control = mclust::emControl(tol = c(1e-12, 1e-12), itmax = c(1e6, 1e6)))
  expect_reference(as.numeric(logLik(candidate)), continued$loglik, reference = "mclust",
                   covers = c("enumerate_lpa()", "get_candidate()"), quantity = "log likelihood, three VVI profiles, continued",
                   tolerance = 1e-8, measure = "relative")
  expect_reference(as.numeric(logLik(candidate_fit(enumeration_fit, n_profiles = 3, model = "VVI"))), continued$loglik,
                   reference = "mclust", covers = "candidate_fit()",
                   quantity = "log likelihood, three VVI profiles, candidate_fit()", tolerance = 1e-8, measure = "relative")
})

test_that("differences, best candidates and the knee follow from mclust's criteria", {
  skip_without_reference("mclust")
  grid <- subset(mclust_grid(), model == "EEI")
  # The decrease in BIC from k - 1 to k profiles within one structure (man/get_deltas.Rd).
  deltas <- data.frame(candidate = grid$candidate, value = grid$bic_individual,
                       delta_profile = c(NA, -diff(grid$bic_individual)))
  expect_reference(subset(get_deltas(enumeration_fit), model == "EEI"), deltas, reference = "mclust",
                   covers = "get_deltas()", quantity = "BIC and its change per added profile, EEI",
                   tolerance = 1e-8, measure = "relative", key = "candidate")
  # Kneedle (Satopaa, Albrecht, Irwin and Raghavan 2011): both axes scaled to [0, 1]; the knee is the
  # point furthest from the chord between the first and the last point.
  x <- (grid$n_profiles - min(grid$n_profiles)) / diff(range(grid$n_profiles))
  y <- (grid$bic_individual - min(grid$bic_individual)) / diff(range(grid$bic_individual))
  last <- length(x)
  distance <- abs((y[last] - y[1]) * x - (x[last] - x[1]) * y + x[last] * y[1] - y[last] * x[1]) /
    sqrt((y[last] - y[1])^2 + (x[last] - x[1])^2)
  knee <- data.frame(line = "1 group class, EEI", knee_candidate = grid$candidate[which.max(distance)],
                     distance = max(distance))
  expect_reference(subset(get_knees(enumeration_fit), direction == "profiles" & line == "1 group class, EEI"), knee,
                   reference = "mclust", covers = "get_knees()", quantity = "knee of the BIC line, EEI",
                   tolerance = 1e-8, measure = "relative", key = "line",
                   note = "The BIC values are mclust's; the knee is computed from them by the kneedle definition.")
  # The best candidate per criterion over the whole grid, the criteria from mclust's likelihoods.
  whole <- mclust_grid()
  whole$log_likelihood[whole$candidate == vvi_three] <- NA
  deviance <- -2 * whole$log_likelihood
  p <- whole$n_parameters
  n <- enumeration_n
  criteria <- list(aic = deviance + 2 * p, kic = deviance + 3 * (p + 1), bic = deviance + p * log(n),
                   sabic = deviance + p * log((n + 2) / 24), caic = deviance + p * (log(n) + 1))
  best <- data.frame(criterion = names(criteria),
                     candidate = vapply(criteria, \(value) whole$candidate[which.min(value)], character(1)),
                     value = vapply(criteria, min, numeric(1), na.rm = TRUE))
  # One row per criterion: the individual-level convention where a criterion has two.
  chosen <- subset(get_best(enumeration_fit), criterion %in% names(criteria) & (is.na(convention) | convention == "individuals"))
  expect_reference(chosen, best, reference = "mclust", covers = "get_best()",
                   quantity = "best candidate and its value per criterion", tolerance = 1e-8, measure = "relative",
                   key = "criterion",
                   note = "Three VVI profiles are left out of mclust's side (a lower local maximum there); latents' best is checked to lie elsewhere.")
  expect_reference(subset(get_results(summary(enumeration_fit), "best"), criterion %in% names(criteria) & (is.na(convention) | convention == "individuals")),
                   best, reference = "mclust", covers = c("summary.multilpa_enumeration()", "get_results.summary_multilpa_enumeration()"),
                   quantity = "best candidate per criterion, summary()", tolerance = 1e-8, measure = "relative", key = "criterion")
})
