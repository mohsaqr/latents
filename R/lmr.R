# The Vuong-Lo-Mendell-Rubin test of K against K + 1 classes (Vuong 1989;
# Lo, Mendell and Rubin 2001), in the form Latent GOLD implements (Vermunt
# 2024). The likelihood-ratio statistic LR = 2 (l1 - l0) is referred to a
# weighted sum of independent chi-square(1) variables whose weights are the
# eigenvalues of
#
#   W = [ -B1 A1^-1   -B10 A0^-1 ]
#       [ B10' A1^-1   B0 A0^-1  ],
#
# with A the Hessian of each model's log likelihood, B the cross-product of
# its per-unit scores and B10 the cross-product of the larger model's scores
# with the smaller one's. The p-value is computed by Imhof's (1961) method.

#' The scores and Hessian of one fitted model, for Vuong's test
#'
#' Read from [parameter_inference()] with `vcov_type = "robust"`, which keeps
#' the per-unit (per-group) scores; its observed information is computed once
#' per session for a fit and reused.
#' @param fit A fitted profile or class model.
#' @return A list of `scores` (units by free coordinates), `information`
#'   (minus the Hessian, in the same coordinates), `log_likelihood` and
#'   `n_parameters`.
#' @noRd
.latents_vuong_parts <- function(fit) {
  inference <- parameter_inference(fit, vcov_type = "robust")
  list(scores = attr(inference, "group_scores"),
       information = attr(inference, "hessian"),
       log_likelihood = fit$log_likelihood, n_parameters = fit$n_parameters)
}

#' The weights of the reference distribution of Vuong's statistic
#'
#' With `A = -information`, `-B A^-1 = B information^-1`, so in terms of the
#' information matrices
#' `W = [ B1 I1^-1, B10 I0^-1 ; -B10' I1^-1, -B0 I0^-1 ]`.
#' For two regular nested models whose information equals the variance of
#' the score, the nonzero eigenvalues are `p1 - p0` ones, and the reference
#' distribution is the familiar chi-square on `p1 - p0` degrees of freedom.
#' @param alternative,null From `.latents_vuong_parts()`, on the same units.
#' @return The real parts of the eigenvalues of `W`.
#' @noRd
.latents_vuong_weights <- function(alternative, null) {
  stopifnot("the two models must be scored on the same units" =
              nrow(alternative$scores) == nrow(null$scores))
  cross_alternative <- crossprod(alternative$scores)
  cross_null <- crossprod(null$scores)
  cross_both <- crossprod(alternative$scores, null$scores)
  inverse_alternative <- solve(alternative$information)
  inverse_null <- solve(null$information)
  weights_matrix <- rbind(
    cbind(cross_alternative %*% inverse_alternative, cross_both %*% inverse_null),
    cbind(-t(cross_both) %*% inverse_alternative, -cross_null %*% inverse_null))
  Re(eigen(weights_matrix, only.values = TRUE)$values)
}

#' Upper tail of a weighted sum of chi-square(1) variables (Imhof 1961)
#'
#' `P(sum_j lambda_j X_j > q)` for independent `X_j ~ chi-square(1)`:
#' `1/2 + (1/pi) int_0^Inf sin(theta(u)) / (u rho(u)) du`, with
#' `theta(u) = sum_j atan(lambda_j u) / 2 - q u / 2` and
#' `rho(u) = prod_j (1 + lambda_j^2 u^2)^(1/4)`. Weights that are zero to
#' rounding contribute nothing and are dropped.
#'
#' The integrand oscillates with half-period `2 pi / q` and, with few
#' weights, decays slowly, which defeats adaptive quadrature over an
#' infinite range. It is integrated half-period by half-period, where it is
#' smooth, and the alternating partial sums are accelerated by Wynn's
#' epsilon algorithm, as QUADPACK's QAWF does for Fourier integrals. The
#' value is accepted when estimates from 40 and 60 more half-periods agree
#' to `1e-9`; otherwise `latents_inaccurate_probability` is raised. The
#' probability is accurate to about `1e-9` in absolute terms: a p-value
#' below that is reported as a value of that order, not its exact size.
#' @param q The statistic.
#' @param lambda The weights.
#' @return A probability.
#' @noRd
.latents_imhof_upper <- function(q, lambda) {
  lambda <- lambda[abs(lambda) > 1e-10 * max(abs(lambda))]
  if (length(lambda) == 0L) return(NA_real_)
  integrand <- function(u) {
    theta <- 0.5 * vapply(u, function(value) sum(atan(lambda * value)), numeric(1)) -
      0.5 * q * u
    rho <- vapply(u, function(value) prod((1 + lambda^2 * value^2)^0.25), numeric(1))
    ifelse(u == 0, 0.5 * (sum(lambda) - q), sin(theta) / (u * rho))
  }
  # Half a period of sin(q u / 2); a statistic near zero oscillates slowly, so
  # the step is capped at the scale of the weights.
  half_period <- min(2 * pi / max(abs(q), 1e-8), 2 * pi / max(abs(lambda)))
  piece <- function(k) {
    stats::integrate(integrand, k * half_period, (k + 1) * half_period,
                     rel.tol = 1e-10, abs.tol = 1e-13)$value
  }
  estimate <- function(pieces) .latents_wynn(cumsum(pieces))
  pieces <- vapply(seq_len(60L) - 1L, piece, numeric(1))
  previous <- estimate(pieces)
  # Each round adds 40 half-periods until two successive estimates agree;
  # the number of rounds is bounded, and reaching the bound is reported.
  for (round in seq_len(25L)) {
    pieces <- c(pieces, vapply(length(pieces) + seq_len(40L) - 1L, piece, numeric(1)))
    current <- estimate(pieces)
    if (abs(current - previous) < 1e-9) {
      # The integral's own rounding can stray outside [0, 1] far in a tail.
      return(min(1, max(0, 0.5 + current / pi)))
    }
    previous <- current
  }
  stop(errorCondition(sprintf(paste(
    "The tail probability of the weighted chi-square distribution did not",
    "settle (statistic %g, %d weights); the p-value cannot be reported."),
    q, length(lambda)), class = "latents_inaccurate_probability", call = NULL))
}

#' The limit of a sequence by Wynn's epsilon algorithm
#'
#' Builds the epsilon table column by column from the partial sums `s`
#' (`e[-1] = 0`, `e[0] = s`, `e[j + 1][k] = e[j - 1][k + 1] + 1 /
#' (e[j][k + 1] - e[j][k])`) and returns the last entry of the last even
#' column, the accelerated limit. Two equal neighbours mean the sequence has
#' already settled there.
#' @param s Partial sums.
#' @return The estimated limit.
#' @noRd
.latents_wynn <- function(s) {
  table <- list(rep(0, length(s) + 1L), s)
  limit <- s[length(s)]
  # The table is built column by column, each from the two before it.
  while (length(table[[length(table)]]) > 2L) {
    current <- table[[length(table)]]
    before <- table[[length(table) - 1L]]
    gaps <- diff(current)
    if (any(gaps == 0)) break
    following <- before[-c(1L, length(before))][seq_along(gaps)] + 1 / gaps
    table <- c(table, list(following))
    if ((length(table) - 2L) %% 2L == 0L) limit <- following[length(following)]
  }
  limit
}

#' One VLMR and adjusted LMR comparison of two fitted candidates
#'
#' The adjusted statistic divides LR by `1 + 1 / (d log n)`, with `d` the
#' difference in the number of parameters and `n` the number of independent
#' units (Lo, Mendell and Rubin 2001, eq. 15, with their three parameters per
#' added univariate normal component replaced by the models' own difference),
#' and is referred to the same distribution.
#' @param alternative,null From `.latents_vuong_parts()`.
#' @return A one-row data frame: `n_units`, `parameters_added`, `lr`,
#'   `vlmr_p`, `adjusted_lr`, `almr_p`, `reference_mean`, `reference_sd`.
#' @noRd
.latents_vuong_row <- function(alternative, null) {
  lambda <- .latents_vuong_weights(alternative, null)
  lr <- 2 * (alternative$log_likelihood - null$log_likelihood)
  added <- alternative$n_parameters - null$n_parameters
  n_units <- nrow(alternative$scores)
  adjusted <- lr / (1 + 1 / (added * log(n_units)))
  data.frame(n_units = n_units, parameters_added = added, lr = lr,
             vlmr_p = .latents_imhof_upper(lr, lambda), adjusted_lr = adjusted,
             almr_p = .latents_imhof_upper(adjusted, lambda),
             reference_mean = sum(lambda), reference_sd = sqrt(2 * sum(lambda^2)))
}

#' The VLMR and adjusted LMR tests between neighbouring candidates of a grid
#'
#' Each converged candidate is compared with the converged candidate of one
#' more profile at the same number of group classes and structure, and with
#' the one of one more group class at the same number of profiles and
#' structure. A comparison that cannot be computed (a candidate on a boundary,
#' a singular information) keeps its row, with `NA` statistics and the reason
#' in `note`; the warnings raised while scoring a candidate are recorded there
#' too.
#' @param x A `multilpa_enumeration`.
#' @return A data frame, one row per comparison: `direction`
#'   (`"profiles"` or `"group_classes"`), `null_candidate`,
#'   `alternative_candidate`, `n_profiles`, `n_group_classes`, `model` (of the
#'   alternative), then the columns of `.latents_vuong_row()` and `note`.
#' @noRd
.latents_lmr_table <- function(x) {
  table <- .latents_named_table(x)
  structures <- ifelse(is.na(table$model), "", table$model)
  usable <- table$converged %in% TRUE
  find <- function(n_profiles, n_group_classes, model) {
    which(usable & table$n_profiles == n_profiles &
            table$n_group_classes == n_group_classes & structures == model)
  }
  pairs <- do.call(rbind, lapply(which(usable), function(at) {
    profile_step <- find(table$n_profiles[at] + 1L, table$n_group_classes[at], structures[at])
    group_step <- find(table$n_profiles[at], table$n_group_classes[at] + 1L, structures[at])
    rbind(if (length(profile_step) == 1L) {
      data.frame(direction = "profiles", null = at, alternative = profile_step)
    }, if (length(group_step) == 1L) {
      data.frame(direction = "group_classes", null = at, alternative = group_step)
    })
  }))
  if (is.null(pairs) || nrow(pairs) == 0L) {
    stop(errorCondition(paste(
      "No two converged candidates of this grid differ by one profile or one",
      "group class, so there is no Lo-Mendell-Rubin test to compute."),
      class = c("latents_no_lmr", "latents_nothing_to_plot"), call = NULL))
  }
  # Each candidate is scored once, whatever the number of comparisons it is in.
  involved <- sort(unique(c(pairs$null, pairs$alternative)))
  scored <- lapply(stats::setNames(involved, involved), function(at) {
    notes <- character()
    parts <- tryCatch(
      withCallingHandlers(.latents_vuong_parts(get_candidate(x, table$candidate[at])),
        warning = function(condition) {
          notes <<- c(notes, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }),
      latents_boundary_fit = function(condition) condition,
      latents_singular_information = function(condition) condition,
      latents_unsupported_inference = function(condition) condition,
      latents_no_converge = function(condition) condition)
    list(parts = parts, notes = notes)
  })
  rows <- lapply(seq_len(nrow(pairs)), function(index) {
    null <- scored[[as.character(pairs$null[index])]]
    alternative <- scored[[as.character(pairs$alternative[index])]]
    failed <- Filter(function(entry) inherits(entry$parts, "condition"),
                     list(null = null, alternative = alternative))
    statistics <- if (length(failed) > 0L) {
      data.frame(n_units = NA_integer_, parameters_added = NA_integer_, lr = NA_real_,
                 vlmr_p = NA_real_, adjusted_lr = NA_real_, almr_p = NA_real_,
                 reference_mean = NA_real_, reference_sd = NA_real_)
    } else .latents_vuong_row(alternative$parts, null$parts)
    reasons <- c(vapply(names(failed), function(role) {
      sprintf("%s model: %s", role, conditionMessage(failed[[role]]$parts))
    }, character(1)), unique(c(null$notes, alternative$notes)))
    alternative_at <- pairs$alternative[index]
    data.frame(direction = pairs$direction[index],
               null_candidate = table$candidate[pairs$null[index]],
               alternative_candidate = table$candidate[alternative_at],
               n_profiles = table$n_profiles[alternative_at],
               n_group_classes = table$n_group_classes[alternative_at],
               model = table$model[alternative_at], statistics,
               note = if (length(reasons) == 0L) NA_character_ else
                 paste(reasons, collapse = " | "),
               stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
