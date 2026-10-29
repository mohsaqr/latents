# Helpers for equivalence/test-features-vlmr.R: reading the kit and the two
# programs' listings, fitting latents' enumeration of a case, and handing a
# fit's estimates to the independent reference (reference.R) as plain numbers.

vlmr_dir <- file.path("features", "vlmr")

#' Evaluate `expr`, muffling only the package's expected fit warnings and its
#' messages (as tests/testthat's quietly(); every other warning stays visible)
vlmr_quietly <- function(expr) {
  expected <- c("latents_unconverged", "latents_failed_starts", "latents_failed_replicates",
                "latents_boundary", "latents_small_classes", "latents_single_level")
  withCallingHandlers(expr,
    warning = function(w) if (inherits(w, expected)) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage"))
}

vlmr_cases <- function() readRDS(file.path(vlmr_dir, "kit", "cases.rds"))

#' Latent GOLD listings are Latin-1 with CRLF line ends
vlmr_lg_lines <- function(name) {
  lines <- readLines(file.path(vlmr_dir, "returned", paste0(name, ".lst")), warn = FALSE)
  sub("\r$", "", iconv(lines, "latin1", "UTF-8"))
}

#' One row per model of a Latent GOLD enumeration listing: `k`, `ll`, `npar`,
#' and for every model after the first the VLMR block: `h0_ll`, `mean`, `sd`,
#' `lr`, `vlmr_p`, `adjusted_lr`, `almr_p` (NA where not printed)
vlmr_read_lg <- function(name) {
  lines <- vlmr_lg_lines(name)
  starts <- grep("^Model[0-9]+ - ", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  rows <- lapply(seq_along(starts), function(index) {
    block <- lines[starts[index]:ends[index]]
    field <- function(label, position = 1L) {
      line <- block[startsWith(block, label)][1L]
      if (is.na(line)) return(NA_real_)
      values <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      as.numeric(values[1L + position])
    }
    k <- as.integer(sub("-Cluster.*$", "", block[2L]))
    data.frame(k = k, ll = field("Log-likelihood (LL)"),
               npar = field("Number of parameters (Npar)"),
               h0_ll = field("LL of H0 model"), mean = field("Mean of VLMR Distr."),
               sd = field("StdDev of VLMR Distr."), lr = field("-2LL Diff (VLMR)"),
               vlmr_p = field("-2LL Diff (VLMR)", 2L),
               adjusted_lr = field("-2LL Diff Adjusted (aVLMR)"),
               almr_p = field("-2LL Diff Adjusted (aVLMR)", 2L))
  })
  do.call(rbind, rows)
}

#' TECH11 of one Mplus run of K classes (K - 1 against K)
vlmr_read_mplus <- function(name, k) {
  path <- file.path(vlmr_dir, "returned", sprintf("%s_k%d.out", name, k))
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  at <- grep("VUONG-LO-MENDELL-RUBIN LIKELIHOOD RATIO TEST", lines)[1L]
  last <- function(line) {
    as.numeric(utils::tail(strsplit(trimws(line), "[[:space:]]+")[[1L]], 1L))
  }
  value <- function(pattern, from = 1L) {
    hit <- grep(pattern, lines)
    hit <- hit[hit >= from][1L]
    if (is.na(hit)) NA_real_ else last(lines[hit])
  }
  h1 <- value("^ +H0 Value", grep("^MODEL FIT INFORMATION", lines)[1L] %||% 1L)
  if (is.na(at)) {
    return(data.frame(k = k, terminated = any(grepl("ESTIMATION TERMINATED NORMALLY", lines)),
                      h1_ll = h1, h0_ll = NA_real_, lr = NA_real_, df = NA_real_,
                      mean = NA_real_, sd = NA_real_, vlmr_p = NA_real_,
                      adjusted_lr = NA_real_, almr_p = NA_real_))
  }
  adjusted_at <- grep("LO-MENDELL-RUBIN ADJUSTED LRT TEST", lines)[1L]
  data.frame(k = k, terminated = any(grepl("ESTIMATION TERMINATED NORMALLY", lines)),
             h1_ll = h1, h0_ll = value("H0 Loglikelihood Value", at),
             lr = value("2 Times the Loglikelihood Difference", at),
             df = value("Difference in the Number of Parameters", at),
             mean = value("^ +Mean ", at), sd = value("^ +Standard Deviation", at),
             vlmr_p = value("^ +P-Value", at),
             adjusted_lr = value("^ +Value ", adjusted_at),
             almr_p = value("^ +P-Value", adjusted_at))
}

#' The arguments that fit a case's enumeration in latents
vlmr_fit_call <- function(case, seed = 1L, n_starts = NULL, data = case$data,
                          range = case$spec$range) {
  spec <- case$spec
  n_starts <- n_starts %||% if (spec$n >= 3000L) 30L else 60L
  common <- list(n_starts = n_starts, seed = seed, tol = 1e-12, max_iter = 5000)
  if (!is.null(spec$weights)) common$weights <- spec$weights
  if (isTRUE(spec$missing)) common$missing <- "fiml"
  if (identical(spec$model, "LCA")) {
    return(vlmr_quietly(do.call(enumerate_lca, c(list(data, spec$vars, n_classes = range), common))))
  }
  vlmr_quietly(do.call(enumerate_lpa, c(list(data, spec$vars, n_profiles = range,
                                        model = spec$model), common)))
}

#' The candidates' log likelihoods of a grid, by number of classes
vlmr_grid_ll <- function(grid) {
  table <- as.data.frame(grid)
  stats::setNames(table$log_likelihood, table$n_profiles)
}

#' Whether each candidate of a grid is bound-active, by number of classes
vlmr_grid_boundary <- function(grid) {
  table <- as.data.frame(grid)
  stats::setNames(table$boundary %in% TRUE, table$n_profiles)
}

#' Fit a case, trying further seeds (with more starts) until every candidate
#' reaches at least the reference maximum less `slack`; the grid reaching the
#' most candidates is kept, with the seed that produced it.
vlmr_fit_case <- function(case, target_ll, slack = 1e-3, seeds = 1:3) {
  attempt <- function(index) {
    grid <- vlmr_fit_call(case, seed = seeds[index],
                          n_starts = if (index == 1L) NULL else 150L)
    ll <- vlmr_grid_ll(grid)[as.character(names(target_ll))]
    list(grid = grid, seed = seeds[index], reached = sum(ll >= target_ll - slack, na.rm = TRUE))
  }
  best <- attempt(1L)
  index <- 1L
  while (best$reached < length(target_ll) && index < length(seeds)) {
    index <- index + 1L
    candidate <- attempt(index)
    if (candidate$reached > best$reached) best <- candidate
  }
  best
}

#' A latents fit's estimates as plain numbers, for the reference
vlmr_estimates <- function(fit, model) {
  k <- fit$n_profiles
  proportions <- get_results(fit, "profile_probabilities")
  proportions <- proportions$probability[order(proportions$profile)]
  if (identical(model, "LCA")) {
    responses <- get_results(fit, "responses")
    second <- subset(responses, as.character(category) == "2")
    second <- second[order(match(second$indicator, fit$vars), second$profile), ]
    return(list(model = model, k = k, proportions = proportions,
                probabilities = matrix(second$probability, k, length(fit$vars))))
  }
  estimates <- list(model = model, k = k, proportions = proportions,
                    means = unname(fit$means), variances = unname(fit$variances))
  if (identical(model, "EEE")) {
    covariances <- subset(get_results(fit, "covariances"), profile == 1)
    vars <- fit$vars
    estimates$covariance <- matrix(covariances$covariance[order(
      match(covariances$indicator_2, vars), match(covariances$indicator, vars))],
      length(vars), length(vars))
  }
  estimates
}

#' The reference comparison of two neighbouring candidates of a latents grid
vlmr_reference_row <- function(grid, case, null_k, draws = 1e6, convention = "sampling") {
  spec <- case$spec
  null_fit <- get_candidate(grid, null_k)
  alternative_fit <- get_candidate(grid, null_k + 1L)
  x <- as.matrix(case$data[spec$vars])
  storage.mode(x) <- "double"
  weights <- if (is.null(spec$weights)) NULL else case$data[[spec$weights]]
  n_units <- if (is.null(weights)) sum(rowSums(!is.na(x)) > 0L) else sum(weights > 0)
  ref_vlmr(vlmr_estimates(alternative_fit, spec$model), vlmr_estimates(null_fit, spec$model),
           x, weights, n_units, draws, convention)
}

#' Latents' own reference weights for one neighbouring pair of a grid
vlmr_latents_weights <- function(grid, null_k) {
  .latents_vuong_weights(.latents_vuong_parts(get_candidate(grid, null_k + 1L)),
                         .latents_vuong_parts(get_candidate(grid, null_k)))
}

#' Fit one case (or read its cached grid) and compute every comparison
#'
#' `VLMR_FIT_CACHE`, when set to a directory, keeps each fitted grid there
#' between runs; it is unset by default, so the equivalence run fits live.
#' Returns `rows`, one per neighbouring pair: latents' table, Latent GOLD's
#' and Mplus's printed values, the reference, and both models' log
#' likelihoods in each program; and `lambda`, latents' and the reference's
#' weights per pair.
vlmr_compute_case <- function(name, case, draws = 1e6) {
  cache <- Sys.getenv("VLMR_FIT_CACHE")
  cached <- if (nzchar(cache)) file.path(cache, paste0(name, ".rds")) else ""
  lg <- if (isTRUE(case$external)) vlmr_read_lg(name) else NULL
  fitted <- if (nzchar(cached) && file.exists(cached)) readRDS(cached) else {
    target <- if (is.null(lg)) NULL else stats::setNames(lg$ll, lg$k)
    result <- if (is.null(target)) list(grid = vlmr_fit_call(case), seed = 1L) else
      vlmr_fit_case(case, target)
    if (nzchar(cached)) saveRDS(result, cached)
    result
  }
  grid <- fitted$grid
  table <- get_results(grid, "lmr")
  ll <- vlmr_grid_ll(grid)
  boundary <- vlmr_grid_boundary(grid)
  rows <- lapply(seq_len(nrow(table)), function(index) {
    row <- table[index, ]
    null_k <- as.integer(row$n_profiles - 1L)
    out <- data.frame(case = name, config = case$config, model = case$spec$model,
                      n = case$spec$n, true_k = case$spec$true_k, null_k = null_k,
                      seed_used = fitted$seed, ll0 = ll[[as.character(null_k)]],
                      ll1 = ll[[as.character(null_k + 1L)]],
                      boundary0 = boundary[[as.character(null_k)]],
                      boundary1 = boundary[[as.character(null_k + 1L)]], row[c(
                        "n_units", "parameters_added", "lr", "vlmr_p", "adjusted_lr", "almr_p",
                        "reference_mean", "reference_sd", "note")])
    if (!is.null(lg)) {
      at <- match(null_k + 1L, lg$k)
      out <- cbind(out, data.frame(lg_ll0 = lg$ll[match(null_k, lg$k)], lg_ll1 = lg$ll[at],
                                   lg_lr = lg$lr[at], lg_vlmr_p = lg$vlmr_p[at],
                                   lg_adjusted_lr = lg$adjusted_lr[at],
                                   lg_almr_p = lg$almr_p[at], lg_mean = lg$mean[at],
                                   lg_sd = lg$sd[at]))
      mplus <- vlmr_read_mplus(name, null_k + 1L)
      if (!is.null(mplus)) {
        names(mplus) <- paste0("mp_", names(mplus))
        out <- cbind(out, mplus[setdiff(names(mplus), "mp_k")])
      }
    }
    list(row = out, null_k = null_k, computed = !is.na(row$lr))
  })
  references <- lapply(rows, function(entry) {
    if (!entry$computed) return(NULL)
    reference <- vlmr_reference_row(grid, case, entry$null_k, draws)
    reference$parts <- NULL
    reference$latents_lambda <- vlmr_latents_weights(grid, entry$null_k)
    if (!is.null(case$spec$weights)) {
      frequency <- vlmr_reference_row(grid, case, entry$null_k, 1e4, "frequency")
      reference$frequency <- c(mean = frequency$reference_mean, sd = frequency$reference_sd,
                               p = frequency$vlmr[["davies"]], almr_p = frequency$almr[["davies"]])
    }
    reference
  })
  # For Latent GOLD: the reference at Latent GOLD's own printed estimates (and
  # with its frequency-weight convention), which measures how far estimates
  # agreeing to 4 decimals move the weights.
  printed <- lapply(rows, function(entry) {
    if (is.null(lg)) return(NULL)
    if (identical(case$spec$model, "LCA")) {
      # Categorical listings are not parsed: only the rounding movement.
      return(if (entry$computed) vlmr_rounding_movement(grid, case, entry$null_k))
    }
    at_lg <- vlmr_reference_at_lg(name, case, entry$null_k)
    if (entry$computed) cbind(at_lg, vlmr_rounding_movement(grid, case, entry$null_k)) else at_lg
  })
  rows0 <- rows
  rows <- Map(function(entry, reference, at_lg) {
    row <- entry$row
    if (!is.null(at_lg)) row <- cbind(row, at_lg)
    if (is.null(reference)) return(row)
    lambda <- sort(reference$lambda, decreasing = TRUE)
    ours <- sort(reference$latents_lambda, decreasing = TRUE)
    cbind(row, data.frame(
      ref_lr = reference$lr, ref_adjusted_lr = reference$adjusted_lr,
      ref_parameters_added = reference$parameters_added, ref_n_units = reference$n_units,
      ref_mean = reference$reference_mean, ref_sd = reference$reference_sd,
      ref_imhof_p = reference$vlmr[["imhof"]], ref_davies_p = reference$vlmr[["davies"]],
      ref_farebrother_p = reference$vlmr[["farebrother"]],
      ref_mc_p = reference$vlmr[["monte_carlo"]],
      ref_imhof_almr_p = reference$almr[["imhof"]], ref_davies_almr_p = reference$almr[["davies"]],
      ref_mc_almr_p = reference$almr[["monte_carlo"]], ref_draws = reference$draws,
      n_weights = length(lambda), n_latents_weights = length(ours),
      lambda_max = max(abs(lambda)),
      lambda_gap = if (length(lambda) == length(ours)) max(abs(lambda - ours)) else NA_real_,
      score_sum = max(reference$score_sums),
      ref_freq_mean = reference$frequency[["mean"]] %||% NA_real_,
      ref_freq_sd = reference$frequency[["sd"]] %||% NA_real_,
      ref_freq_p = reference$frequency[["p"]] %||% NA_real_,
      ref_freq_almr_p = reference$frequency[["almr_p"]] %||% NA_real_))
  }, rows, references, printed)
  lambda <- lapply(Filter(Negate(is.null), references), function(reference) {
    list(latents = reference$latents_lambda, reference = reference$lambda)
  })
  names(lambda) <- sprintf("%s:%d", name, vapply(Filter(function(entry) entry$computed, rows0),
                                                  `[[`, integer(1), "null_k"))
  list(rows = vlmr_bind(rows), seed = fitted$seed, lambda = lambda)
}

#' Map fixed-seed cases without forking an unsafe BLAS runtime
#'
#' Fitting and reference Monte Carlo each set explicit seeds. Serial execution
#' therefore preserves the case streams and can also be requested through
#' `options(latents.vlmr.serial = TRUE)`.
#' @param cases The case names or inputs to map.
#' @param fun A function of one case.
#' @param cores Maximum number of workers on a fork-safe platform.
#' @return A list in input order.
vlmr_map_cases <- function(cases, fun, cores = 2L) {
  if (isTRUE(getOption("latents.vlmr.serial", FALSE)) ||
      identical(.Platform$OS.type, "windows") || !latents:::.latents_fork_safe()) {
    lapply(cases, fun)
  } else {
    parallel::mclapply(cases, fun, mc.cores = cores, mc.preschedule = FALSE)
  }
}

#' Every case, on at most `cores` safe workers; results bound into one table
vlmr_compute_all <- function(names, cores = 2L, draws = 1e6) {
  cases <- vlmr_cases()
  results <- vlmr_map_cases(names, function(name) {
    tryCatch(vlmr_compute_case(name, cases[[name]], draws),
             error = function(condition) list(error = conditionMessage(condition)))
  }, cores = cores)
  names(results) <- names
  failed <- Filter(function(result) !is.null(result$error), results)
  succeeded <- Filter(function(result) is.null(result$error), results)
  list(table = vlmr_bind(lapply(succeeded, `[[`, "rows")),
       lambda = do.call(c, unname(lapply(succeeded, `[[`, "lambda"))), errors = vapply(failed, `[[`, character(1), "error"))
}

#' Bind data frames on the union of their columns (missing ones NA)
vlmr_bind <- function(frames) {
  columns <- unique(unlist(lapply(frames, names)))
  table <- do.call(rbind, lapply(frames, function(frame) {
    frame[setdiff(columns, names(frame))] <- NA
    frame[columns]
  }))
  rownames(table) <- NULL
  table
}

#' P(a X_m + b Y_k > q) for X ~ chi-square(m), Y ~ chi-square(k), a > 0, by
#' one-dimensional quadrature of the convolution (y = t^2 removes the
#' chi-square(1) singularity at zero); exact to the quadrature's 1e-13
vlmr_two_group_tail <- function(q, a, m, b, k) {
  density <- function(t) if (k == 1) 2 * stats::dnorm(t) else 2 * t * stats::dchisq(t^2, k)
  inner <- function(t) {
    stats::pchisq(pmax((q - b * t^2) / a, 0), m, lower.tail = FALSE) * density(t)
  }
  # Y's tail beyond t = 40 sqrt(k) + 40 carries no representable mass; with
  # b > 0 the integrand is zero past the kink sqrt(q / b), where P(Y > q / b)
  # takes over. Breakpoints keep each segment within the density's scale.
  upper <- 40 * sqrt(k) + 40
  tail <- 0
  if (b > 0) {
    kink <- sqrt(max(q, 0) / b)
    tail <- stats::pchisq(kink^2, k, lower.tail = FALSE)
    upper <- min(upper, kink)
  }
  if (upper <= 0) return(tail)
  breaks <- sort(unique(c(0, pmin(sqrt(k) * c(0.5, 1, 2, 4, 8), upper), upper)))
  tail + sum(vapply(seq_len(length(breaks) - 1L), function(index) {
    stats::integrate(inner, breaks[index], breaks[index + 1L], rel.tol = 1e-13,
                     abs.tol = 1e-16, subdivisions = 5000L)$value
  }, numeric(1)))
}

#' The stress set for .latents_imhof_upper(): one row per weight vector and
#' statistic, with an exact probability where one exists and CompQuadForm's
#' Davies (accuracy 1e-10) otherwise
vlmr_imhof_stress <- function(seed = 20261020L) {
  set.seed(seed)
  row <- function(family, q, lambda, exact, exact_kind = "exact") {
    latents <- .latents_imhof_upper(q, lambda)
    data.frame(family = family, n_weights = length(lambda), q = q,
               weight_range = max(abs(lambda)) / min(abs(lambda)), latents = latents,
               reference = exact, reference_kind = exact_kind,
               abs_error = abs(latents - exact),
               # Relative errors are reported where the probability exceeds 1e-6.
               rel_error = if (exact > 1e-6) abs(latents - exact) / exact else NA_real_)
  }
  one <- lapply(seq_len(60), function(index) {
    scale <- 10^stats::runif(1L, -6, 6)
    q <- scale * 10^stats::runif(1L, -4, 1.9)
    row("one weight", q, scale, stats::pchisq(q / scale, 1, lower.tail = FALSE))
  })
  negative <- lapply(seq_len(20), function(index) {
    scale <- 10^stats::runif(1L, -3, 3)
    q <- -scale * 10^stats::runif(1L, -3, 1.5)
    row("one negative weight", q, -scale, stats::pchisq(-q / scale, 1))
  })
  equal <- lapply(seq_len(60), function(index) {
    scale <- 10^stats::runif(1L, -6, 6)
    df <- sample(2:30, 1L)
    q <- scale * stats::qchisq(stats::runif(1L, 1e-6, 1 - 1e-6), df)
    row("equal weights", q, rep(scale, df), stats::pchisq(q / scale, df, lower.tail = FALSE))
  })
  two <- lapply(seq_len(60), function(index) {
    a <- 10^stats::runif(1L, -2, 2)
    b <- sample(c(-1, 1), 1L) * 10^stats::runif(1L, -2, 2)
    q <- abs(a + b) * stats::rexp(1L) + stats::runif(1L, 0, 3) * a
    row("two weights", q, c(a, b), vlmr_two_group_tail(q, a, 1, b, 1))
  })
  f_ratio <- lapply(seq_len(60), function(index) {
    m <- sample(1:10, 1L)
    k <- sample(1:10, 1L)
    ratio <- 10^stats::runif(1L, -4, 4)
    scale <- 10^stats::runif(1L, -6, 6)
    row("mixed signs, q = 0", 0, scale * c(rep(1, m), rep(-ratio, k)),
        stats::pf(ratio * k / m, m, k, lower.tail = FALSE))
  })
  grouped <- lapply(seq_len(40), function(index) {
    m <- sample(5:45, 1L)
    a <- 10^stats::runif(1L, -1, 1)
    b <- sample(c(-1, 1), 1L) * 10^stats::runif(1L, -1, 1)
    q <- max(0, (m * a + (50 - m) * b) + stats::rnorm(1L, 0, 3) * sqrt(m * a^2 + (50 - m) * b^2))
    row("50 weights in two groups", q, c(rep(a, m), rep(b, 50 - m)),
        vlmr_two_group_tail(q, a, m, b, 50 - m))
  })
  many <- lapply(seq_len(40), function(index) {
    lambda <- stats::rnorm(50) * 10^stats::runif(50, -1, 1)
    q <- max(0, sum(lambda) + stats::rnorm(1L, 0, 1) * sqrt(2 * sum(lambda^2)))
    row("50 random weights", q, lambda,
        CompQuadForm::davies(q, lambda, acc = 1e-11, lim = 1e7)$Qq, "davies")
  })
  dynamic <- lapply(seq_len(30), function(index) {
    m <- sample(1:5, 1L)
    k <- sample(1:5, 1L)
    small <- 10^stats::runif(1L, -8, -6)
    large <- 10^stats::runif(1L, 3, 4)
    q <- large * stats::qchisq(stats::runif(1L, 0.01, 0.99), k)
    row("dynamic range 1e-8 to 1e4", q, c(rep(small, m), rep(large, k)),
        vlmr_two_group_tail(q, large, k, small, m))
  })
  extreme_q <- lapply(seq_len(30), function(index) {
    df <- sample(1:6, 1L)
    q <- if (index %% 2L == 0L) 10^stats::runif(1L, -12, -6) else 10^stats::runif(1L, 2, 3)
    row(if (index %% 2L == 0L) "tiny q" else "huge q", q, rep(1, df),
        stats::pchisq(q, df, lower.tail = FALSE))
  })
  scaled <- lapply(seq_len(60), function(index) {
    lambda <- stats::rnorm(6) * 10^stats::runif(6, -1, 1)
    q <- max(0, sum(lambda) + stats::rnorm(1L) * sqrt(2 * sum(lambda^2)))
    size <- c(1e-6, 1e6)[1L + index %% 2L]
    row(sprintf("joint scaling by %g", size), q * size, lambda * size,
        .latents_imhof_upper(q, lambda), "unscaled latents")
  })
  do.call(rbind, c(one, negative, equal, two, f_ratio, grouped, many, dynamic, extreme_q,
                   scaled))
}

#' Latent GOLD's printed estimates of one continuous model of a listing, as
#' reference.R takes them (`NULL` for a categorical model)
#'
#' With `parameters=first`, class j's mean is the intercept plus its effect
#' and the shares are the softmax of the printed logits. A model is
#' degenerate when a printed share is below 1 / n or a variance prints 0.
vlmr_lg_estimates <- function(name, k, model, vars, n) {
  if (identical(model, "LCA")) return(NULL)
  lines <- vlmr_lg_lines(name)
  starts <- grep("^Model[0-9]+ - ", lines)
  at <- which(as.integer(sub("-Cluster.*$", "", lines[starts + 1L])) == k)
  block <- lines[starts[at]:(c(starts[-1L] - 1L, length(lines))[at])]
  first <- grep("^Cluster\\(1\\)\t<-\t1\t", block)[1L]
  fields <- lapply(strsplit(block[first:length(block)], "\t", fixed = TRUE),
                   function(field) field[nzchar(field)])
  value <- function(test) {
    hit <- Filter(test, fields)
    if (length(hit) == 0L) NA_real_ else as.numeric(hit[[1L]][length(hit[[1L]]) - 0L])
  }
  coefficient <- function(pattern) {
    hit <- Filter(function(field) identical(field[seq_along(pattern)], pattern), fields)[1L]
    as.numeric(hit[[1L]][length(pattern) + 1L])
  }
  logits <- vapply(seq_len(k), function(j) coefficient(c(sprintf("Cluster(%d)", j), "<-", "1")),
                   numeric(1))
  shares <- exp(logits - max(logits))
  shares <- shares / sum(shares)
  means <- vapply(vars, function(v) {
    coefficient(c(v, "<-", "1")) +
      vapply(seq_len(k), function(j) coefficient(c(v, "<-", sprintf("Cluster(%d)", j))),
             numeric(1))
  }, numeric(k))
  means <- matrix(means, k, length(vars))
  variances <- if (identical(model, "VVI")) {
    matrix(vapply(vars, function(v) vapply(seq_len(k), function(j) {
      coefficient(c(v, "|", sprintf("Cluster(%d)", j)))
    }, numeric(1)), numeric(k)), k, length(vars))
  } else {
    matrix(vapply(vars, function(v) {
      hit <- Filter(function(field) identical(field[1L], v) && length(field) >= 2L &&
                      !field[2L] %in% c("<-", "|", "<->", "Continuous"), fields)[1L]
      as.numeric(hit[[1L]][2L])
    }, numeric(1)), k, length(vars), byrow = TRUE)
  }
  estimates <- list(model = model, k = k, proportions = shares, means = means,
                    variances = variances,
                    degenerate = any(shares < 1 / n) || any(variances <= 0))
  if (identical(model, "EEE")) {
    covariance <- diag(variances[1L, ], length(vars))
    pairs <- utils::combn(length(vars), 2L)
    covariance[t(pairs)] <- covariance[t(pairs[2:1, ])] <- apply(pairs, 2L, function(pair) {
      coefficient(c(vars[pair[1L]], "<->", vars[pair[2L]]))
    })
    estimates$covariance <- covariance
  }
  estimates
}

#' The reference evaluated at Latent GOLD's printed estimates of one pair
#'
#' Returns one row: whether either printed model is degenerate (a share below
#' 1 / n or a zero variance), the log likelihoods there, and the moments and
#' p-values of the reference (frequency weights, as Latent GOLD's caseweight).
vlmr_reference_at_lg <- function(name, case, null_k) {
  spec <- case$spec
  alternative <- vlmr_lg_estimates(name, null_k + 1L, spec$model, spec$vars, spec$n)
  null <- vlmr_lg_estimates(name, null_k, spec$model, spec$vars, spec$n)
  degenerate <- isTRUE(alternative$degenerate) || isTRUE(null$degenerate)
  empty <- data.frame(lg_degenerate = degenerate, lgpt_ll0 = NA_real_, lgpt_ll1 = NA_real_,
                      lgpt_mean = NA_real_, lgpt_sd = NA_real_, lgpt_p = NA_real_,
                      lgpt_almr_p = NA_real_)
  if (degenerate) return(empty)
  x <- as.matrix(case$data[spec$vars])
  storage.mode(x) <- "double"
  weights <- if (is.null(spec$weights)) NULL else case$data[[spec$weights]]
  # Latent GOLD counts every row it reads (`missing includeall`).
  reference <- ref_vlmr(alternative, null, x, weights, n_units = nrow(x), draws = 1e4,
                        convention = "frequency")
  a <- reference$parts$alternative
  b <- reference$parts$null
  data.frame(lg_degenerate = FALSE, lgpt_ll0 = b$log_likelihood, lgpt_ll1 = a$log_likelihood,
             lgpt_mean = reference$reference_mean, lgpt_sd = reference$reference_sd,
             lgpt_p = reference$vlmr[["davies"]], lgpt_almr_p = reference$almr[["davies"]])
}

#' How far the reference moves when latents' estimates move by Latent GOLD's
#' printed rounding: six draws of uniform +-5e-5 on every mean, variance
#' (covariance), response probability and class logit, the largest change of each quantity. Uses
#' Latent GOLD's frequency-weight convention and its row count n.
vlmr_rounding_movement <- function(grid, case, null_k, draws = 6L) {
  spec <- case$spec
  x <- as.matrix(case$data[spec$vars])
  storage.mode(x) <- "double"
  weights <- if (is.null(spec$weights)) NULL else case$data[[spec$weights]]
  estimates <- list(alternative = vlmr_estimates(get_candidate(grid, null_k + 1L), spec$model),
                    null = vlmr_estimates(get_candidate(grid, null_k), spec$model))
  summary <- function(alternative, null) {
    reference <- ref_vlmr(alternative, null, x, weights, n_units = nrow(x), draws = 1e3,
                          convention = "frequency")
    c(mean = reference$reference_mean, sd = reference$reference_sd,
      p = reference$vlmr[["davies"]], almr_p = reference$almr[["davies"]])
  }
  base <- do.call(summary, estimates)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
    assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(case$seed + null_k)
  jitter <- function(values) values + stats::runif(length(values), -5e-5, 5e-5)
  perturb <- function(estimate) {
    logits <- jitter(log(estimate$proportions))
    estimate$proportions <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
    if (identical(estimate$model, "LCA")) {
      estimate$probabilities[] <- pmin(pmax(jitter(estimate$probabilities), 1e-12), 1 - 1e-12)
      return(estimate)
    }
    estimate$means[] <- jitter(estimate$means)
    if (identical(estimate$model, "EEI")) {
      estimate$variances[] <- matrix(jitter(estimate$variances[1L, ]), estimate$k,
                                     ncol(estimate$variances), byrow = TRUE)
    } else estimate$variances[] <- jitter(estimate$variances)
    if (!is.null(estimate$covariance)) {
      upper <- jitter(estimate$covariance[upper.tri(estimate$covariance, diag = TRUE)])
      estimate$covariance[upper.tri(estimate$covariance, diag = TRUE)] <- upper
      estimate$covariance[lower.tri(estimate$covariance)] <-
        t(estimate$covariance)[lower.tri(estimate$covariance)]
    }
    estimate
  }
  moved <- vapply(seq_len(draws), function(index) {
    abs(summary(perturb(estimates$alternative), perturb(estimates$null)) - base)
  }, numeric(4))
  data.frame(round_mean = max(moved[1L, ]), round_sd = max(moved[2L, ]),
             round_p = max(moved[3L, ]), round_almr_p = max(moved[4L, ]))
}
