# Functions for compare.R: growth mixtures in mixture_regression() against
# genuine Mplus 9 Demo runs kept in JStats (tests/fixtures/mplus/; copied to
# returned/ with their data): latent class growth (growth-factor variances
# fixed at zero) at one to three classes on two datasets, and growth mixtures
# with a class-specific random intercept and slope covariance at one to three
# classes. Every class has its own residual variance, equal across waves
# (Mplus labels (r1), (r2), ...), which is latents' default structure.
#
# Tolerances, fixed before the comparison: Mplus prints three decimals, so
# half a unit (5e-4) plus 1e-6 of the value; counts exactly.

# delta_terms() of the wider Latent GOLD kit, found from the project root or
# from equivalence/.
source(if (dir.exists(file.path("equivalence", "latentgold-wider"))) {
  file.path("equivalence", "latentgold-wider", "compare-functions.R")
} else file.path("latentgold-wider", "compare-functions.R"))

.mplus_cases <- function() {
  list(
    list(name = "lcga_ml_k1", folder = "lcga-ml", file = "lcga-k1-ml.out", data = "lcga.dat", classes = 1L, random = FALSE),
    list(name = "lcga_ml_k2", folder = "lcga-ml", file = "lcga-k2-ml.out", data = "lcga.dat", classes = 2L, random = FALSE),
    list(name = "lcga_ml_k3", folder = "lcga-ml", file = "lcga-k3-ml.out", data = "lcga.dat", classes = 3L, random = FALSE),
    list(name = "runner_lcga_k1", folder = "runner-lcga", file = "lcga-k1.out", data = "lcga.dat", classes = 1L, random = FALSE, estimator = "MLR"),
    list(name = "runner_lcga_k2", folder = "runner-lcga", file = "lcga-k2.out", data = "lcga.dat", classes = 2L, random = FALSE, estimator = "MLR"),
    list(name = "runner_lcga_k3", folder = "runner-lcga", file = "lcga-k3.out", data = "lcga.dat", classes = 3L, random = FALSE, estimator = "MLR"),
    list(name = "gmm_ml_k2", folder = "gmm-ml", file = "gmm-k2-ml.out", data = "gmm.dat", classes = 2L, random = TRUE),
    list(name = "runner_gmm_k1", folder = "runner-gmm-growth", file = "gmm-k1.out", data = "gmm.dat", classes = 1L, random = TRUE, estimator = "MLR"),
    list(name = "runner_gmm_k2", folder = "runner-gmm-growth", file = "gmm-k2.out", data = "gmm.dat", classes = 2L, random = TRUE, estimator = "MLR"),
    list(name = "runner_gmm_k3", folder = "runner-gmm-growth", file = "gmm-k3.out", data = "gmm.dat", classes = 3L, random = TRUE, estimator = "MLR"))
}

#' Long data from an Mplus wide file, waves at time 0, 1, ...
#' @noRd
.mplus_long <- function(path) {
  wide <- as.matrix(utils::read.table(path))
  data.frame(id = rep(seq_len(nrow(wide)), each = ncol(wide)),
             time = rep(seq_len(ncol(wide)) - 1L, times = nrow(wide)), y = as.vector(t(wide)))
}

#' Fit statistics and the per-class growth parameters an Mplus output prints
#' @return A list: `log_likelihood`, `n_parameters`, `bic`, and `parameters`
#'   (data frame: `term`, `estimate`, `standard_error`, printed text kept).
read_mplus_growth <- function(path) {
  lines <- readLines(path, warn = FALSE)
  value <- function(pattern) {
    hit <- grep(pattern, lines, value = TRUE)[1L]
    utils::tail(strsplit(trimws(hit), "[[:space:]]+")[[1L]], 1L)
  }
  start <- grep("^MODEL RESULTS", lines)[1L]
  end <- grep("^(QUALITY OF NUMERICAL RESULTS|RESULTS IN PROBABILITY SCALE|TECHNICAL 1 OUTPUT|LOGISTIC REGRESSION ODDS RATIO RESULTS)",
              lines)
  end <- end[end > start][1L]
  block <- lines[seq.int(start, end - 1L)]
  class <- 1L
  section <- ""
  rows <- list()
  invisible(lapply(block, function(line) {
    if (grepl("^Latent Class [0-9]+", line)) {
      class <<- as.integer(sub("^Latent Class ([0-9]+).*", "\\1", line))
      section <<- ""
      return(NULL)
    }
    if (grepl("^Categorical Latent Variables", line)) {
      class <<- 0L
      section <<- ""
      return(NULL)
    }
    if (grepl("^ [A-Z]", line)) {
      section <<- trimws(gsub("[[:space:]]+", " ", line))
      return(NULL)
    }
    fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    if (length(fields) == 5L && grepl("^[A-Z]", fields[[1L]]) && grepl("^-?[0-9.]+$", fields[[2L]])) {
      rows[[length(rows) + 1L]] <<- data.frame(
        term = if (class == 0L) sprintf("%s %s", section, fields[[1L]]) else
          sprintf("class %d %s %s", class, section, fields[[1L]]),
        estimate = fields[[2L]], standard_error = fields[[3L]], stringsAsFactors = FALSE)
    }
  }))
  parameters <- do.call(rbind, rows)
  # Fixed loadings and intercepts print 999.000 as their test statistic; only
  # estimated terms are kept.
  parameters <- parameters[as.numeric(parameters$standard_error) > 0 &
                             !grepl("\\| ", parameters$term) &
                             !grepl("Intercepts", parameters$term), , drop = FALSE]
  parameters <- parameters[!duplicated(sub("(Residual Variances) Y[0-9]+", "\\1", parameters$term)) |
                             !grepl("Residual Variances", parameters$term), , drop = FALSE]
  list(log_likelihood = value("^ +H0 Value"), n_parameters = as.integer(value("^Number of Free Parameters")),
       bic = value("^ +Bayesian \\(BIC\\)"), parameters = parameters)
}

fit_mplus_case <- function(case, data) {
  arguments <- list(y ~ time, data, case$classes, id = "id", class_level = "group",
                    n_starts = if (case$classes == 1L) 1L else 50L, seed = 1L, tol = 1e-12,
                    max_iter = 100000L)
  if (case$random) arguments <- c(arguments, list(random = ~time, random_covariance = "varying"))
  if (identical(case$estimator, "MLR")) arguments$vcov_type <- "robust"
  do.call(mixture_regression, arguments)
}

#' Mplus's terms from a latents growth fit under one class order
#' @noRd
mplus_growth_map <- function(perm, random, waves_first) {
  function(theta) {
    k <- length(perm)
    gamma <- function(j) if (j == 1L) 0 else theta[[sprintf("membership.class_%d.(Intercept)", j)]]
    by_class <- lapply(seq_len(k), function(c) {
      j <- perm[[c]]
      out <- c(theta[[sprintf("coefficient.class_%d.(Intercept)", j)]],
               theta[[sprintf("coefficient.class_%d.time", j)]],
               exp(2 * theta[[sprintf("log_sigma.class_%d", j)]]))
      names(out) <- sprintf("class %d %s", c, c("Means I", "Means S",
                                                 sprintf("Residual Variances Y%d", waves_first)))
      if (random) {
        a <- theta[[sprintf("random.log_cholesky.class_%d.(Intercept):(Intercept)", j)]]
        b <- theta[[sprintf("random.log_cholesky.class_%d.time:time", j)]]
        lower <- theta[[sprintf("random.cholesky.class_%d.time:(Intercept)", j)]]
        out <- c(out, stats::setNames(c(exp(2 * a), lower^2 + exp(2 * b), exp(a) * lower),
                                      sprintf("class %d %s", c, c("Variances I", "Variances S",
                                                                  "S WITH I"))))
      }
      out
    })
    # Mplus's class logits are against its last class.
    logits <- if (k > 1L) stats::setNames(vapply(seq_len(k - 1L), function(c)
      gamma(perm[[c]]) - gamma(perm[[k]]), numeric(1)), sprintf("Means C#%d", seq_len(k - 1L)))
    c(unlist(by_class), logits)
  }
}

.printed_rows <- function(case, quantities, printed, obtained, kind = "commercial program") {
  digits <- nchar(sub("^-?[0-9]*[.]?", "", printed))
  value <- as.numeric(printed)
  tolerance <- 0.5 * 10^(-digits) + 1e-6 * abs(value)
  data.frame(case = case, reference = "Mplus 9 Demo (JStats fixtures)", quantity = quantities,
             reference_value = value, latents = obtained, difference = abs(obtained - value),
             tolerance = tolerance, n_values = 1L, kind = kind,
             agrees = abs(obtained - value) <= tolerance)
}

compare_mplus_case <- function(case, directory) {
  folder <- file.path(directory, "returned", case$folder)
  data <- .mplus_long(file.path(folder, case$data))
  listing <- read_mplus_growth(file.path(folder, case$file))
  fit <- fit_mplus_case(case, data)
  theta <- stats::coef(fit)
  first <- grep("Residual Variances", listing$parameters$term, value = TRUE)[1L]
  wave <- as.integer(sub(".*Y([0-9]+)$", "\\1", first))
  perms <- .permutations(case$classes)
  maps <- lapply(perms, mplus_growth_map, random = case$random, waves_first = wave)
  printed_means <- grep("Means I$", listing$parameters$term)
  gaps <- vapply(maps, function(map) {
    value <- map(theta)
    max(abs(value[listing$parameters$term[printed_means]] -
              as.numeric(listing$parameters$estimate[printed_means])))
  }, numeric(1))
  # The runner series uses Mplus's default MLR: the person-level sandwich
  # without a small-sample factor. latents' robust covariance carries the CR1
  # factor G / (G - 1) (G persons), removed here.
  covariance <- stats::vcov(fit)
  if (identical(case$estimator, "MLR")) {
    persons <- length(unique(data$id))
    covariance <- covariance * (persons - 1) / persons
  }
  terms <- delta_terms(maps[[which.min(gaps)]], theta, covariance)
  at <- match(listing$parameters$term, terms$term)
  stopifnot("every Mplus term has a latents counterpart" = !anyNA(at))
  rbind(
    .printed_rows(case$name, c("log likelihood", "BIC"), c(listing$log_likelihood, listing$bic),
                  c(as.numeric(stats::logLik(fit)), stats::BIC(fit))),
    data.frame(case = case$name, reference = "Mplus 9 Demo (JStats fixtures)",
               quantity = "number of parameters", reference_value = listing$n_parameters,
               latents = fit$n_parameters, difference = abs(fit$n_parameters - listing$n_parameters),
               tolerance = 0, n_values = 1L, kind = "commercial program",
               agrees = fit$n_parameters == listing$n_parameters),
    .printed_rows(case$name, paste("estimate", listing$parameters$term), listing$parameters$estimate,
                  terms$estimate[at]),
    .printed_rows(case$name, paste("s.e.", listing$parameters$term), listing$parameters$standard_error,
                  terms$standard_error[at]))
}

.permutations <- function(k) {
  if (k == 1L) return(list(1L))
  unlist(lapply(seq_len(k), function(first) lapply(.permutations(k - 1L), function(rest)
    c(first, setdiff(seq_len(k), first)[rest]))), recursive = FALSE)
}

compare_all_mplus_growth <- function(directory) {
  out <- do.call(rbind, lapply(.mplus_cases(), compare_mplus_case, directory = directory))
  out$fn <- "growth"
  out
}
