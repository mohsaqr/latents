# Functions for compare.R: lpa() and lca() against genuine Stata 19 output
# printed in its manual (Stata [SEM] Example 52g, latent profile models;
# Example 50g, a latent class model) and against the published model-fit table
# of Masyn (2013), Table 25.12, on the same diabetes data. Masyn's adjusted
# LMR p-values (Mplus TECH11) are not compared: TECH11 fits the K-class null
# by deleting the first class of the K + 1 solution, a different test from
# latents' (and Latent GOLD's) Vuong statistic between the two fitted models.
#
# Tolerances, fixed before the comparison: half a unit of the last printed
# digit plus 1e-6 of the value (both programs' convergence and numerical
# derivatives; Stata prints about seven significant digits).

source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

.printed_row <- function(case, source, quantity, printed, obtained, kind) {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  value <- as.numeric(printed)
  tolerance <- 0.5 * 10^(-digits) + 1e-6 * abs(value)
  data.frame(case = case, reference = source, quantity = quantity, reference_value = value,
             latents = obtained, difference = abs(obtained - value), tolerance = tolerance,
             n_values = 1L, kind = kind, agrees = abs(obtained - value) <= tolerance)
}

#' A printed Stata table: estimates and standard errors against latents terms
#' @param printed Named list: term -> c(estimate, s.e.) as printed text.
#' @param terms delta_terms() output under the chosen class matching.
#' @noRd
.table_rows <- function(case, source, printed, terms, kind) {
  at <- match(names(printed), terms$term)
  stopifnot("every printed term has a latents counterpart" = !anyNA(at))
  do.call(rbind, lapply(seq_along(printed), function(i) rbind(
    .printed_row(case, source, paste("estimate", names(printed)[[i]]), printed[[i]][[1L]],
                 terms$estimate[[at[[i]]]], kind),
    .printed_row(case, source, paste("s.e.", names(printed)[[i]]), printed[[i]][[2L]],
                 terms$standard_error[[at[[i]]]], kind))))
}

#' Stata's terms of a latent profile fit under one class order
#' @noRd
stata_lpa_map <- function(perm, vars, shared) {
  function(theta) {
    k <- length(perm)
    share <- function(j) theta[[sprintf("profile.probability.profile_%d.group_class_1", j)]]
    mean <- function(j, v) theta[[sprintf("measurement.mean.profile_%d.%s", j, v)]]
    covariance <- function(j, a, b) {
      owner <- if (shared) "shared" else sprintf("profile_%d", j)
      if (shared) theta[[sprintf("measurement.variance.shared.%s", a)]] else
        theta[[sprintf("measurement.covariance.%s.%s:%s", owner, b, a)]]
    }
    pairs <- if (shared) NULL else utils::combn(vars, 2L, simplify = FALSE)
    c(stats::setNames(vapply(seq_len(k)[-1L], function(c) log(share(perm[c]) / share(perm[1L])),
                             numeric(1)), sprintf("logit %d.C", seq_len(k)[-1L])),
      unlist(lapply(seq_len(k), function(c) c(
        stats::setNames(vapply(vars, function(v) mean(perm[c], v), numeric(1)),
                        sprintf("class %d mean %s", c, vars)),
        stats::setNames(vapply(vars, function(v) covariance(perm[c], v, v), numeric(1)),
                        sprintf("class %d var %s", c, vars)),
        if (!shared) stats::setNames(vapply(pairs, function(p) covariance(perm[c], p[[1L]], p[[2L]]),
                                            numeric(1)),
                                     sprintf("class %d cov %s %s", c, vapply(pairs, `[[`, "", 1L),
                                             vapply(pairs, `[[`, "", 2L)))))))
  }
}

#' The class order whose means are closest to Stata's printed means
#' @noRd
.closest_perm <- function(maps, theta, printed) {
  mean_terms <- grep(" mean ", names(printed), value = TRUE)
  target <- vapply(printed[mean_terms], function(x) as.numeric(x[[1L]]), numeric(1))
  gaps <- vapply(maps, function(map) max(abs(map(theta)[mean_terms] - target)), numeric(1))
  which.min(gaps)
}

fit_reaven <- function(data, classes, model) {
  lpa(data, c("glucose", "insulin", "sspg"), n_profiles = classes, model = model, n_starts = 50L,
      seed = 1L, tol = 1e-15, max_iter = 200000L)
}

compare_stata_lpa <- function(data) {
  vars <- c("glucose", "insulin", "sspg")
  source <- "Stata 19 [SEM] Example 52g (gsem)"
  eei <- list(`logit 2.C` = c("-1.541025", ".2205682"),
              `class 1 mean glucose` = c("41.22237", "1.298051"),
              `class 1 mean insulin` = c("20.98005", "1.000974"),
              `class 1 mean sspg` = c("14.96579", ".6868081"),
              `class 1 var glucose` = c("191.5596", "23.83815"),
              `class 1 var insulin` = c("119.0542", "14.00336"),
              `class 1 var sspg` = c("55.91283", "6.713667"),
              `class 2 mean glucose` = c("115.7123", "2.849914"),
              `class 2 mean insulin` = c("7.553144", "2.160949"),
              `class 2 mean sspg` = c("34.5529", "1.53117"))
  vvv <- list(`logit 2.C` = c("-.8853513", ".2386536"), `logit 3.C` = c("-.612664", ".2260018"),
              `class 1 mean glucose` = c("35.68584", ".5741752"),
              `class 1 mean insulin` = c("16.58066", ".6204724"),
              `class 1 mean sspg` = c("10.49755", ".5833606"),
              `class 1 var glucose` = c("19.30952", "3.932547"),
              `class 1 var insulin` = c("26.7354", "4.494093"),
              `class 1 var sspg` = c("18.71079", "3.970509"),
              `class 1 cov glucose insulin` = c("3.456027", "2.942391"),
              `class 1 cov glucose sspg` = c("5.474303", "2.811729"),
              `class 1 cov insulin sspg` = c("7.995803", "3.020304"),
              `class 2 mean glucose` = c("47.66176", "1.492718"),
              `class 2 mean insulin` = c("34.35203", "3.00337"),
              `class 2 mean sspg` = c("24.414", ".7395383"),
              `class 2 var glucose` = c("53.21326", "15.56547"),
              `class 2 var insulin` = c("228.6332", "59.03553"),
              `class 2 var sspg` = c("13.75515", "3.838523"),
              `class 2 cov glucose insulin` = c("40.02875", "23.12762"),
              `class 2 cov glucose sspg` = c(".7294854", "5.48065"),
              `class 2 cov insulin sspg` = c("-5.743169", "11.4943"),
              `class 3 mean glucose` = c("93.92473", "6.985336"),
              `class 3 mean insulin` = c("10.37614", "1.123135"),
              `class 3 mean sspg` = c("28.4787", "1.94975"),
              `class 3 var glucose` = c("1279.011", "312.6774"),
              `class 3 var insulin` = c("36.38521", "9.26287"),
              `class 3 var sspg` = c("113.3239", "27.67628"),
              `class 3 cov glucose insulin` = c("-163.4383", "47.637"),
              `class 3 cov glucose sspg` = c("276.9206", "81.60543"),
              `class 3 cov insulin sspg` = c("-25.4313", "11.66564"))
  run <- function(case, classes, model, printed, loglik) {
    fit <- fit_reaven(data, classes, model)
    theta <- stats::coef(fit)
    maps <- lapply(permutations(classes), stata_lpa_map, vars = vars, shared = identical(model, "EEI"))
    terms <- delta_terms(maps[[.closest_perm(maps, theta, printed)]], theta, stats::vcov(fit))
    rbind(.printed_row(case, source, "log likelihood", loglik, as.numeric(stats::logLik(fit)),
                       "commercial program"),
          .table_rows(case, source, printed, terms, "commercial program"))
  }
  rbind(run("reaven_EEI_k2", 2L, "EEI", eei, "-1702.5542"),
        run("reaven_VVV_k3", 3L, "VVV", vvv, "-1536.6409"))
}

#' Masyn (2013) Table 25.12: log likelihood, parameters, BIC by structure
#' @noRd
compare_masyn <- function(data) {
  source <- "Masyn (2013), Oxford Handbook of Quantitative Methods, Table 25.12"
  table <- data.frame(
    model = rep(c("EEI", "VVI", "EEE", "VVV"), c(5L, 5L, 4L, 3L)),
    classes = c(1:5, 1:5, 1:4, 1:3),
    loglik = c("-1820.68", "-1702.55", "-1653.24", "-1606.30", "-1578.21",
               "-1820.68", "-1641.95", "-1562.48", "-1544.10", "-1528.73",
               "-1730.40", "-1666.63", "-1628.86", "-1591.84",
               "-1730.40", "-1590.57", "-1536.64"),
    npar = c(6, 10, 14, 18, 22, 6, 13, 20, 27, 34, 9, 13, 17, 21, 9, 19, 29),
    bic = c("3671.22", "3454.88", "3376.15", "3302.18", "3265.90",
            "3671.22", "3348.60", "3224.49", "3222.57", "3226.67",
            "3505.60", "3397.95", "3342.33", "3288.19",
            "3505.60", "3275.69", "3217.61"), stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(table)), function(i) {
    row <- table[i, ]
    fit <- fit_reaven(data, row$classes, row$model)
    case <- sprintf("reaven_%s_k%d", row$model, row$classes)
    rbind(.printed_row(case, source, "log likelihood", row$loglik, as.numeric(stats::logLik(fit)),
                       "published"),
          data.frame(case = case, reference = source, quantity = "number of parameters",
                     reference_value = row$npar, latents = fit$n_parameters,
                     difference = abs(fit$n_parameters - row$npar), tolerance = 0, n_values = 1L,
                     kind = "published", agrees = fit$n_parameters == row$npar),
          .printed_row(case, source, "BIC", row$bic, stats::BIC(fit), "published"))
  }))
}

#' Stata [SEM] Example 50g: two latent classes of the Stouffer-Toby items
#' @noRd
compare_stata_lca <- function(data) {
  source <- "Stata 19 [SEM] Example 50g (gsem)"
  items <- c("accident", "play", "insurance", "stock")
  fit <- lca(data, items, n_classes = 2L, n_starts = 50L, seed = 1L, tol = 1e-15,
             max_iter = 200000L)
  theta <- stats::coef(fit)
  map <- function(perm) function(theta) {
    p <- function(j, v) theta[[sprintf("measurement.response.profile_%d.%s:1", j, v)]]
    share <- function(j) theta[[sprintf("profile.probability.profile_%d.group_class_1", j)]]
    c(`logit 2.C` = log(share(perm[2L]) / share(perm[1L])),
      unlist(lapply(1:2, function(c) stats::setNames(
        vapply(items, function(v) stats::qlogis(p(perm[c], v)), numeric(1)),
        sprintf("class %d logit %s", c, items)))),
      `class 1 probability` = share(perm[1L]), `class 2 probability` = share(perm[2L]))
  }
  printed <- list(`logit 2.C` = c("-.9482041", ".2886333"),
                  `class 1 logit accident` = c(".9128742", ".1974695"),
                  `class 1 logit play` = c("-.7099072", ".2249096"),
                  `class 1 logit insurance` = c("-.6014307", ".2123096"),
                  `class 1 logit stock` = c("-1.880142", ".3337665"),
                  `class 2 logit accident` = c("4.983017", "3.745987"),
                  `class 2 logit play` = c("2.747366", "1.165853"),
                  `class 2 logit insurance` = c("2.534582", ".9644841"),
                  `class 2 logit stock` = c("1.203416", ".5361735"),
                  `class 1 probability` = c(".7207539", ".0580926"),
                  `class 2 probability` = c(".2792461", ".0580926"))
  perms <- list(c(1L, 2L), c(2L, 1L))
  gaps <- vapply(perms, function(perm) abs(map(perm)(theta)[["class 1 probability"]] - 0.7207539),
                 numeric(1))
  terms <- delta_terms(map(perms[[which.min(gaps)]]), theta, stats::vcov(fit))
  rbind(.printed_row("stouffer_toby_k2", source, "log likelihood", "-504.46767",
                     as.numeric(stats::logLik(fit)), "commercial program"),
        .table_rows("stouffer_toby_k2", source, printed, terms, "commercial program"))
}

compare_all_stata <- function(directory) {
  reaven <- utils::read.csv(file.path(directory, "reaven-miller.csv"))
  stouffer <- utils::read.csv(file.path(directory, "stouffer-toby.csv"))
  lpa_rows <- rbind(compare_stata_lpa(reaven), compare_masyn(reaven))
  lpa_rows$fn <- "lpa"
  lca_rows <- compare_stata_lca(stouffer)
  lca_rows$fn <- "lca"
  rbind(lpa_rows, lca_rows)
}
