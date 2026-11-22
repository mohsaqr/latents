# The edition of DESCRIPTION (Config/testthat/edition: 3), stated here because a re-run
# of this file outside the package finds no DESCRIPTION and would fall back to edition 2.
testthat::local_edition(3)

# get_tna() and get_group_tna() of multilpa() fits made with `time` against tna
# itself: the networks are tna's relative transition model of the fit's modal
# profile sequences (one row per group, one column per occasion), and
# get_group_tna() splits the sequences by their modal group class. The modal
# sequences are taken by Bayes' rule from coef() in base R (as in
# test-accessors-transitions.R) and handed to tna::tna() and tna::group_tna().
# Ratios of integer counts: 1e-12.

transition_panel <- local({
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, globalenv()), add = TRUE)
  set.seed(5)
  n <- 150L
  occasions <- 4L
  kind <- 1L + (stats::runif(n) < 0.5)
  z <- matrix(stats::rnorm(n * occasions), n)
  states <- do.call(cbind, Reduce(function(previous, t) {
    move <- stats::runif(n) < stats::plogis(c(-1.8, -0.2)[kind] + 0.8 * z[, t])
    ifelse(move, 3L - previous, previous)
  }, seq_len(occasions)[-1L], init = 1L + (stats::runif(n) < 0.4), accumulate = TRUE))
  s <- as.vector(t(states))
  data.frame(id = rep(seq_len(n), each = occasions), time = rep(seq_len(occasions), n),
             z = as.vector(t(z)), w = rep(stats::rnorm(n), each = occasions),
             y1 = stats::rnorm(length(s), c(-1, 1.2)[s], c(0.8, 1)[s]),
             y2 = stats::rnorm(length(s), c(0, 1.5)[s], 0.9))
})
panel_n <- 150L
panel_y1 <- matrix(transition_panel$y1, panel_n, 4L, byrow = TRUE)
panel_y2 <- matrix(transition_panel$y2, panel_n, 4L, byrow = TRUE)
panel_w <- matrix(transition_panel$w, panel_n, 4L, byrow = TRUE)[, 1L]
panel_states <- c("profile_1", "profile_2")

rows_log_sum_exp <- function(values) {
  top <- apply(values, 1L, max)
  top + log(rowSums(exp(values - top)))
}

# Modal profile of every row (sequences x occasions) and modal class of every
# sequence, by Bayes' rule: `log_share[[h]]` the log profile probabilities in
# group class h, `log_class` the log class probabilities (sequences x classes).
modal_sequences <- function(fit, log_share, log_class) {
  estimates <- coef(fit)
  at <- \(name, ...) estimates[[sprintf(name, ...)]]
  log_density <- lapply(1:2, \(k) {
    stats::dnorm(panel_y1, at("measurement.mean.profile_%d.y1", k), sqrt(at("measurement.variance.profile_%d.y1", k)),
                 log = TRUE) +
      stats::dnorm(panel_y2, at("measurement.mean.profile_%d.y2", k), sqrt(at("measurement.variance.profile_%d.y2", k)),
                   log = TRUE)
  })
  per_class <- lapply(1:2, \(h) {
    terms <- lapply(1:2, \(k) log_share[[h]][k] + log_density[[k]])
    row_total <- log(exp(terms[[1L]]) + exp(terms[[2L]]))
    list(terms = terms, row_total = row_total, group = log_class[, h] + rowSums(row_total))
  })
  group_total <- rows_log_sum_exp(cbind(per_class[[1L]]$group, per_class[[2L]]$group))
  class_posterior <- lapply(per_class, \(part) exp(part$group - group_total))
  profile_posterior <- lapply(1:2, \(k) {
    Reduce(`+`, lapply(1:2, \(h) class_posterior[[h]] * exp(per_class[[h]]$terms[[k]] - per_class[[h]]$row_total)))
  })
  states <- ifelse(profile_posterior[[1L]] >= profile_posterior[[2L]], 1L, 2L)
  sequences <- as.data.frame(matrix(panel_states[states], nrow(states)))
  names(sequences) <- sprintf("T%d", seq_len(ncol(sequences)))
  list(sequences = sequences, class = ifelse(class_posterior[[1L]] >= class_posterior[[2L]], 1L, 2L))
}

network_table <- function(network, class) {
  weights <- network$weights[panel_states, panel_states]
  data.frame(class = class, from = rep(panel_states, 2), to = rep(panel_states, each = 2),
             weight = as.vector(weights))
}

check_networks_with_tna <- function(fit, modal, covers) {
  setting <- class(fit)[1L]
  reference <- tna::tna(modal$sequences)
  expect_reference(network_table(get_tna(fit), 0L), network_table(reference, 0L), reference = "tna",
                   covers = c("get_tna()", covers[["pooled"]]),
                   quantity = sprintf("pooled network weights, %s", setting), tolerance = 1e-12,
                   measure = "relative", key = c("class", "from", "to"))
  expect_reference(get_tna(fit)$inits[panel_states], reference$inits[panel_states], reference = "tna",
                   covers = c("get_tna()", covers[["pooled"]]),
                   quantity = sprintf("pooled initial distribution, %s", setting), tolerance = 1e-12,
                   measure = "relative")
  groups <- get_group_tna(fit)
  classes <- sort(unique(modal$class))
  grouped <- tna::group_tna(modal$sequences, group = sprintf("class_%d", modal$class))
  expect_reference(
    do.call(rbind, lapply(seq_along(classes), \(i) network_table(groups[[i]], classes[i]))),
    do.call(rbind, lapply(classes, \(h) network_table(grouped[[sprintf("class_%d", h)]], h))),
    reference = "tna", covers = c("get_group_tna()", covers[["grouped"]]),
    quantity = sprintf("network weights per modal group class, %s", setting), tolerance = 1e-12,
    measure = "relative", key = c("class", "from", "to"))
}

test_that("get_tna() and get_group_tna() of a multilpa() fit are tna's networks of its modal sequences", {
  skip_without_reference("tna")
  fit <- quietly(multilpa(transition_panel, c("y1", "y2"), "id", n_profiles = 2, n_group_classes = 2,
                          time = "time", n_starts = 2, seed = 1, tol = 1e-10))
  estimates <- coef(fit)
  log_share <- lapply(1:2, \(h) log(estimates[sprintf("profile.probability.profile_%d.group_class_%d", 1:2, h)]))
  log_class <- matrix(log(estimates[sprintf("group.probability.group_class_%d", 1:2)]), panel_n, 2L, byrow = TRUE)
  check_networks_with_tna(fit, modal_sequences(fit, log_share, log_class),
                          c(pooled = "get_tna.multilpa()", grouped = "get_group_tna.multilpa()"))
})

test_that("get_tna() and get_group_tna() of a multilpa() fit with a group covariate are tna's networks", {
  skip_without_reference("tna")
  fit <- quietly(multilpa(transition_panel, c("y1", "y2"), "id", n_profiles = 2, n_group_classes = 2,
                          time = "time", group_covariates = "w", n_starts = 2, seed = 1, tol = 1e-10))
  estimates <- coef(fit)
  log_share <- lapply(1:2, \(h) {
    logits <- c(estimates[[sprintf("profile.coefficient.profile_1.group_class_%d", h)]], 0)
    logits - log(sum(exp(logits)))
  })
  eta <- estimates[["group.coefficient.group_class_1.(Intercept)"]] +
    estimates[["group.coefficient.group_class_1.w"]] * panel_w
  log_class <- cbind(stats::plogis(eta, log.p = TRUE), stats::plogis(-eta, log.p = TRUE))
  check_networks_with_tna(fit, modal_sequences(fit, log_share, log_class),
                          c(pooled = "get_tna.multilpa_covariates()", grouped = "get_group_tna.multilpa_covariates()"))
})
