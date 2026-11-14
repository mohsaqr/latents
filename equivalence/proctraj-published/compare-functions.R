# Functions for compare.R: latent class growth models (mixture_regression()
# with one class per person and a residual standard deviation shared by the
# classes) against SAS PROC TRAJ (Jones, Nagin and Roeder 2001) as printed in
# Andruff, Carraro, Thompson, Gaudreau and Louvet (2009), Latent class growth
# modelling: a tutorial, Tutorials in Quantitative Methods for Psychology
# 5(1), 11-24: disengagement coping of 107 persons at three occasions (time
# 0, 1, 2; Appendix A), PROC TRAJ's censored normal model with MAX 5. No
# score lies at a bound (1.00 to 4.63), so its likelihood is the normal one.
#
# Compared: Figure 3 (one quadratic trajectory) and Figure 4 (three) in full
# -- every estimate and standard error, Sigma, the group percentages, the log
# likelihood L, AIC and both BICs -- and the BICs the text reports for one
# linear and two quadratic trajectories. PROC TRAJ's criteria are on the
# log-likelihood scale: AIC = L - p, BIC = L - p log(N) / 2 with N the
# observations (321) or the persons (107).
#
# Not compared: Figure 5 and Table 2 (quadratic, linear and quadratic
# trajectories: mixture_regression() fits one formula in every class) and the
# four-trajectory BIC (the text does not give its polynomial orders).
#
# Tolerances, fixed before the comparison: half a unit of the last printed
# digit plus 1e-6 of the value.
#
# Standard errors: every error both figures print is latents' observed-
# information error times sqrt(N / (N - q)), N = 321 observations and q the
# parameters other than Sigma (1.004706 with q = 3 in Figure 3, 1.017588 with
# q = 11 in Figure 4), to the printed digits: a degrees-of-freedom factor on
# the residual variance. The factor was measured on these outputs, not taken
# from PROC TRAJ's documentation; latents' errors are scaled by it before the
# comparison.

.proctraj_source <- "SAS PROC TRAJ (Andruff et al. 2009)"

#' One comparison row of a printed value
#' @noRd
.proctraj_row <- function(case, quantity, printed, obtained) {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  value <- as.numeric(printed)
  tolerance <- 0.5 * 10^(-digits) + 1e-6 * abs(value)
  difference <- abs(obtained - value)
  data.frame(case = case, reference = .proctraj_source, quantity = quantity,
             reference_value = value, latents = obtained, difference = difference,
             tolerance = tolerance, n_values = 1L, kind = "commercial program",
             agrees = difference <= tolerance, fn = "growth")
}

#' The disengagement data in long form (one row per person and occasion)
#' @noRd
read_disengagement <- function(directory) {
  wide <- utils::read.csv(file.path(directory, "disen.csv"))
  long <- stats::reshape(wide, direction = "long", varying = c("dise1", "dise2", "dise3"),
                         v.names = "disengagement", timevar = "time", times = 0:2, idvar = "id")
  long <- long[order(long$id, long$time), ]
  rownames(long) <- NULL
  long
}

#' A latent class growth model: one class per person, a shared Sigma
#' @noRd
fit_trajectories <- function(data, k, formula = disengagement ~ time + I(time^2)) {
  mixture_regression(formula, data, n_classes = k, id = "id", class_level = "group",
                     variance = "equal", n_starts = 50L, seed = 1L, tol = 1e-12,
                     max_iter = 10000L, vcov_type = "observed")
}

#' PROC TRAJ's criteria: L, AIC, BIC on observations and on persons
#' @param printed Named character vector with any of `L`, `AIC`, `BIC_obs`,
#'   `BIC_persons`.
#' @noRd
.proctraj_criteria <- function(case, fit, printed) {
  log_likelihood <- as.numeric(stats::logLik(fit))
  p <- attr(stats::logLik(fit), "df")
  obtained <- c(L = log_likelihood, AIC = log_likelihood - p,
                BIC_obs = log_likelihood - p * log(321) / 2,
                BIC_persons = log_likelihood - p * log(107) / 2)
  labels <- c(L = "log likelihood L", AIC = "AIC", BIC_obs = "BIC (N=321)",
              BIC_persons = "BIC (N=107)")
  do.call(rbind, lapply(names(printed), function(name)
    .proctraj_row(case, labels[[name]], printed[[name]], obtained[[name]])))
}

#' A printed PROC TRAJ table: estimates and errors of each group, Sigma and
#' the group percentages
#' @param printed List, one element per printed group: named `c(estimate,
#'   s.e.)` text for Intercept, Linear, Quadratic.
#' @noRd
.proctraj_table_rows <- function(case, fit, printed, sigma, groups = NULL) {
  coefficients <- as.data.frame(get_results(fit, "coefficients"))
  classes <- as.data.frame(get_results(fit, "classes"))
  q <- attr(stats::logLik(fit), "df") - 1L
  factor <- sqrt(321 / (321 - q))
  terms <- c(Intercept = "(Intercept)", Linear = "time", Quadratic = "I(time^2)")
  # PROC TRAJ's group g is the latents class whose intercept is closest.
  intercepts <- coefficients[coefficients$term == "(Intercept)", ]
  class_of <- vapply(printed, function(group) {
    intercepts$class[[which.min(abs(intercepts$estimate - as.numeric(group$Intercept[[1L]])))]]
  }, character(1))
  stopifnot("each printed group is a different class" = !anyDuplicated(class_of))
  rows <- do.call(rbind, lapply(seq_along(printed), function(g) {
    do.call(rbind, lapply(names(printed[[g]]), function(parameter) {
      row <- coefficients[coefficients$class == class_of[[g]] &
                            coefficients$term == terms[[parameter]], ]
      rbind(.proctraj_row(case, sprintf("group %d %s estimate", g, parameter),
                          printed[[g]][[parameter]][[1L]], row$estimate),
            .proctraj_row(case, sprintf("group %d %s s.e.", g, parameter),
                          printed[[g]][[parameter]][[2L]], row$std_error * factor))
    }))
  }))
  sigma_rows <- rbind(.proctraj_row(case, "Sigma estimate", sigma[[1L]], classes$sigma[[1L]]),
                      .proctraj_row(case, "Sigma s.e.", sigma[[2L]],
                                    classes$sigma_std_error[[1L]] * factor))
  group_rows <- if (length(groups)) do.call(rbind, lapply(seq_along(groups), function(g) {
    row <- classes[classes$class == class_of[[g]], ]
    rbind(.proctraj_row(case, sprintf("group %d membership (%%) estimate", g), groups[[g]][[1L]],
                        100 * row$prior),
          .proctraj_row(case, sprintf("group %d membership (%%) s.e.", g), groups[[g]][[2L]],
                        100 * row$prior_std_error * factor))
  }))
  rbind(rows, sigma_rows, group_rows)
}

compare_proctraj <- function(directory) {
  data <- read_disengagement(directory)
  stopifnot("107 persons at three occasions" = nrow(data) == 321L &&
              length(unique(data$id)) == 107L)
  one <- fit_trajectories(data, 1L)
  linear <- fit_trajectories(data, 1L, disengagement ~ time)
  two <- fit_trajectories(data, 2L)
  three <- fit_trajectories(data, 3L)
  rbind(
    # Figure 3: one quadratic trajectory.
    .proctraj_criteria("proctraj_k1", one, c(L = "-300.45", AIC = "-304.45", BIC_obs = "-311.99",
                                             BIC_persons = "-309.79")),
    .proctraj_table_rows("proctraj_k1", one,
                         list(list(Intercept = c("2.08570", "0.05992"),
                                   Linear = c("0.02794", "0.15278"),
                                   Quadratic = c("0.00776", "0.07339"))),
                         sigma = c("0.61695", "0.02446")),
    # The text: one linear trajectory and two quadratic trajectories.
    .proctraj_criteria("proctraj_k1_linear", linear, c(BIC_persons = "-307.79")),
    .proctraj_criteria("proctraj_k2", two, c(BIC_persons = "-289.05")),
    # Figure 4: three quadratic trajectories.
    .proctraj_criteria("proctraj_k3", three, c(L = "-255.59", AIC = "-267.59", BIC_obs = "-290.22",
                                               BIC_persons = "-283.63")),
    .proctraj_table_rows("proctraj_k3", three,
                         list(list(Intercept = c("1.78076", "0.06816"),
                                   Linear = c("-0.33513", "0.17409"),
                                   Quadratic = c("0.16754", "0.08295")),
                              list(Intercept = c("2.30806", "0.07786"),
                                   Linear = c("0.28371", "0.16913"),
                                   Quadratic = c("-0.10033", "0.08264")),
                              list(Intercept = c("3.17624", "0.26002"),
                                   Linear = c("1.46432", "0.64282"),
                                   Quadratic = c("-0.69299", "0.30563"))),
                         sigma = c("0.43853", "0.01937"),
                         groups = list(c("47.45755", "6.99427"), c("49.33060", "6.91892"),
                                       c("3.21185", "1.89869"))))
}
