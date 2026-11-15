# An independent reference for Wald inference in a homogeneous latent
# transition model (`lta()`, class `multilpa_transitions`) under a constrained
# covariance structure. It shares nothing with the package's inference code:
# its own forward recursion over each sequence (written from the definition,
# in log space, with its own Gaussian density), the structure charts of
# helper-structure-reference.R (Helmert log shapes, Givens-angle orientations,
# a trace-zero matrix logarithm), softmax logits against the *first* category
# for every distribution (the package uses the last), and numerical rather
# than analytic derivatives. At a maximum the standard errors of the natural
# parameters do not depend on the chart, so the two must agree.

# Log likelihood of a homogeneous first-order transition model with one or
# more group classes. `sequences` is a list of observation-row indices, one
# per group, in occasion order; `data_x` is the indicator matrix.
reference_lta_loglik <- function(data_x, sequences, means, sigma, initial, transition,
                                 group_probabilities) {
  k <- nrow(means)
  # A row with missing indicators contributes the marginal density of what it
  # observed (missing at random); rows are grouped by their observed pattern.
  observed <- !is.na(data_x)
  patterns <- apply(observed, 1L, function(row) paste(as.integer(row), collapse = ""))
  log_density <- matrix(0, nrow(data_x), k)
  for (pattern in unique(patterns)) {  # one block of rows per pattern
    rows <- which(patterns == pattern)
    seen <- observed[rows[1L], ]
    if (!any(seen)) next
    log_density[rows, ] <- vapply(seq_len(k), function(profile) {
      factor <- chol(sigma[[profile]][seen, seen, drop = FALSE])
      residuals <- sweep(data_x[rows, seen, drop = FALSE], 2L, means[profile, seen], "-")
      solved <- backsolve(factor, t(residuals), transpose = TRUE)
      -0.5 * colSums(solved^2) - sum(log(diag(factor))) - sum(seen) / 2 * log(2 * base::pi)
    }, numeric(length(rows)))
  }
  log_sum <- function(v) {
    top <- max(v)
    top + log(sum(exp(v - top)))
  }
  per_group <- vapply(sequences, function(rows) {
    by_class <- vapply(seq_along(group_probabilities), function(class) {
      log_transition <- log(matrix(transition[, , class], k, k))
      alpha <- log(initial[class, ]) + log_density[rows[1L], ]
      for (row in rows[-1L]) {  # the forward recursion, one occasion at a time
        alpha <- vapply(seq_len(k), function(to) {
          log_sum(alpha + log_transition[, to])
        }, numeric(1)) + log_density[row, ]
      }
      log(group_probabilities[class]) + log_sum(alpha)
    }, numeric(1))
    log_sum(by_class)
  }, numeric(1))
  sum(per_group)
}

# The reference chart of a transition fit: means and the structure's own
# spread chart (from helper-structure-reference.R, through a cross-sectional
# view of the fit), then group-class, initial and transition logits against
# the first category.
reference_lta_chart <- function(fit) {
  k <- fit$n_profiles
  h <- fit$n_group_classes
  shim <- list(covariance_structure = fit$covariance_structure,
               n_profiles = k, n_group_classes = 1L,
               means = unname(fit$means), variances = unname(fit$variances),
               covariances = if (is.null(fit$covariances)) NULL else unname(fit$covariances),
               profile_probabilities = matrix(rep(1 / k, k), 1L),
               group_probabilities = 1)
  spread <- reference_structure_chart(shim)
  n_spread <- length(spread$start) - (k - 1L)
  first_logits <- function(p) log(p[-1L]) - log(p[1L])
  softmax_first <- function(logits) {
    v <- c(0, logits)
    exp(v - max(v)) / sum(exp(v - max(v)))
  }
  pieces <- list(
    measurement = spread$start[seq_len(n_spread)],
    group = first_logits(fit$group_probabilities),
    initial = unlist(lapply(seq_len(h), function(class) {
      first_logits(fit$initial_probabilities[class, ])
    })),
    transition = unlist(lapply(seq_len(h), function(class) {
      unlist(lapply(seq_len(k), function(from) {
        first_logits(fit$transition_probabilities[from, , class])
      }))
    })))
  widths <- lengths(pieces)
  decode <- function(psi) {
    at <- cumsum(c(0L, widths))
    block <- function(i) psi[at[i] + seq_len(widths[i])]
    measurement <- spread$decode(c(block(1L), numeric(k - 1L)))
    initial <- matrix(block(3L), h, k - 1L, byrow = TRUE)
    moves <- block(4L)
    transition <- array(0, c(k, k, h))
    for (class in seq_len(h)) for (from in seq_len(k)) {  # fill one row at a time
      at_row <- ((class - 1L) * k + from - 1L) * (k - 1L)
      transition[from, , class] <- softmax_first(moves[at_row + seq_len(k - 1L)])
    }
    list(means = measurement$means, sigma = measurement$sigma,
         group_probabilities = softmax_first(block(2L)),
         initial = t(apply(initial, 1L, softmax_first)),
         transition = transition)
  }
  list(start = unlist(pieces, use.names = FALSE), decode = decode)
}

# The natural parameters in parameter_inference()'s order: means, the spread
# block (shared once for EII, per profile otherwise), group-class shares,
# initial probabilities by group class, then transition probabilities by group
# class, origin and destination.
reference_lta_natural <- function(fit, parameters) {
  d <- ncol(fit$means)
  k <- fit$n_profiles
  h <- fit$n_group_classes
  ellipsoidal <- !endsWith(fit$covariance_structure, "I")
  rows <- if (identical(fit$covariance_structure, "EII")) 1L else seq_len(k)
  lower <- which(lower.tri(matrix(0, d, d), diag = TRUE))
  spread <- unlist(lapply(rows, function(p) {
    s <- parameters$sigma[[p]]
    if (ellipsoidal) s[lower] else diag(s)
  }))
  c(as.vector(t(parameters$means)), spread, parameters$group_probabilities,
    as.vector(t(parameters$initial)),
    unlist(lapply(seq_len(h), function(class) {
      as.vector(t(matrix(parameters$transition[, , class], k, k)))
    })))
}

# Natural-scale standard errors of a transition fit from the reference chart.
reference_lta_structure_se <- function(fit, data, hessian = reference_hessian,
                                       jacobian = reference_jacobian) {
  chart <- reference_lta_chart(fit)
  ordered <- order(data[[fit$id]], data[[fit$time]])
  data_x <- as.matrix(data[fit$continuous])
  sequences <- unname(split(ordered, data[[fit$id]][ordered]))
  loglik <- function(psi) {
    p <- chart$decode(psi)
    reference_lta_loglik(data_x, sequences, p$means, p$sigma, p$initial, p$transition,
                         p$group_probabilities)
  }
  information <- -hessian(loglik, chart$start)
  natural <- function(psi) reference_lta_natural(fit, chart$decode(psi))
  map <- jacobian(natural, chart$start)
  covariance <- map %*% solve(information, t(map))
  list(standard_error = sqrt(pmax(diag(covariance), 0)),
       estimate = natural(chart$start),
       log_likelihood = loglik(chart$start))
}
