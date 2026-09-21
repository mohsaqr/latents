# Readers for the latent transition fixtures copied from the JStats
# repository (tests/fixtures/ there). Provenance, retrieval date and md5 of
# each file are in fixtures/jstats/PROVENANCE.md.

jstats_fixture <- function(name) {
  path <- equivalence_fixture("jstats", name)
  connection <- if (grepl("[.]gz$", name)) gzfile(path) else path
  jsonlite::fromJSON(connection)
}

# A Gaussian scenario with its univariate case brought to the multivariate
# shapes: JStats writes d = 1 data as subjects x occasions and the means and
# standard deviations as length-K vectors, which jsonlite reads as such.
jstats_gaussian <- function(scenario) {
  k <- scenario$k
  d <- scenario$d
  panel <- scenario$data
  if (length(dim(panel)) == 2L) dim(panel) <- c(dim(panel), 1L)
  scenario$data <- panel
  scenario$mu <- matrix(scenario$mu, k, d)
  scenario$sd <- matrix(scenario$sd, k, d)
  scenario
}

# JStats stores panel data as a subjects x occasions x indicators array. The
# long frame has one row per subject and occasion, indicators named `prefix1`,
# `prefix2`, ... in array order.
jstats_panel <- function(panel, prefix) {
  stopifnot("`panel` must be a three-dimensional array" =
              is.array(panel) && length(dim(panel)) == 3L)
  sizes <- dim(panel)
  vars <- paste0(prefix, seq_len(sizes[3L]))
  columns <- lapply(seq_len(sizes[3L]), \(m) as.vector(t(panel[, , m])))
  names(columns) <- vars
  data.frame(id = rep(seq_len(sizes[1L]), each = sizes[2L]),
             time = rep(seq_len(sizes[2L]), sizes[1L]),
             columns)
}

# A frequency-weighted pattern table (one row per response pattern, the last
# column its count) expanded to one row per subject and occasion. Columns run
# occasion-major: occasion 1 indicators 1..M, occasion 2 indicators 1..M, ...
jstats_patterns <- function(patterns, n_occasions, n_indicators, prefix) {
  stopifnot("`patterns` must have occasions x indicators + 1 columns" =
              ncol(patterns) == n_occasions * n_indicators + 1L)
  counts <- patterns[, ncol(patterns)]
  people <- patterns[rep(seq_len(nrow(patterns)), counts),
                     seq_len(n_occasions * n_indicators), drop = FALSE]
  panel <- array(people, c(nrow(people), n_indicators, n_occasions))
  jstats_panel(aperm(panel, c(1L, 3L, 2L)), prefix)
}

# The profile permutation that maps multilpa's labels onto the reference's:
# element k is the multilpa profile matching reference profile k, found by
# ranking one profile-specific quantity in both solutions.
match_profiles <- function(ours, theirs) {
  stopifnot(length(ours) == length(theirs), !anyDuplicated(ours),
            !anyDuplicated(theirs))
  order(ours)[rank(theirs)]
}

# Marginal log likelihood of a first-order homogeneous hidden Markov model by
# brute-force summation over every latent path. Deliberately shares no code
# with the package's forward-backward recursion; feasible because the fixtures
# have K^T <= 81 paths. `log_emission(y, k)` returns the log density of one
# occasion's indicator vector under profile k.
path_log_likelihood <- function(data, vars, initial, transition, log_emission) {
  k <- length(initial)
  occasions <- sort(unique(data$time))
  paths <- as.matrix(expand.grid(rep(list(seq_len(k)), length(occasions))))
  subjects <- split(data[vars], data$id)
  per_subject <- vapply(subjects, function(frame) {
    y <- as.matrix(frame)
    emission <- vapply(seq_len(k), \(profile) {
      vapply(seq_len(nrow(y)), \(t) log_emission(y[t, ], profile), numeric(1))
    }, numeric(nrow(y)))
    steps <- if (ncol(paths) > 1L) {
      rowSums(matrix(log(transition[cbind(as.vector(paths[, -ncol(paths)]),
                                          as.vector(paths[, -1L]))]),
                     nrow(paths)))
    } else 0
    reached <- rowSums(matrix(emission[cbind(rep(seq_len(ncol(paths)), each = nrow(paths)),
                                             as.vector(paths))], nrow(paths)))
    log_path <- log(initial[paths[, 1L]]) + steps + reached
    top <- max(log_path)
    top + log(sum(exp(log_path - top)))
  }, numeric(1))
  sum(per_subject)
}

gaussian_emission <- function(means, sds) {
  function(y, profile) sum(stats::dnorm(y, means[profile, ], sds[profile, ], log = TRUE))
}

binary_emission <- function(success) {
  function(y, profile) sum(stats::dbinom(y, 1L, success[profile, ], log = TRUE))
}

# multilpa's response probabilities as a profiles x indicators matrix of
# P(indicator = 1).
success_probabilities <- function(fit) {
  vapply(fit$response_probabilities, \(block) block[, "1"],
         numeric(fit$n_profiles))
}
