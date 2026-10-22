# Reading the pair-covariance kit's Latent GOLD and Mplus outputs, and
# fitting each case with latents. Sourced by
# equivalence/test-features-pair-covariance.R (run from equivalence/).
# Every class-level quantity is returned as one tidy table:
# class, quantity (proportion, mean, variance, covariance), term, estimate,
# standard_error, with the reference's own class labels.

pc_root <- function(...) file.path("features", "pair-covariance", ...)
pc_cases <- function() readRDS(pc_root("kit", "cases.rds"))

#' latents fitted as the programs were asked to fit the case; the best of
#' several seeds (each with its own random starts).
pc_fit <- function(case, seeds = 1:3, n_starts = 20L) {
  spec <- case$spec
  fits <- lapply(seeds, function(seed) {
    arguments <- list(case$data, spec$vars, spec$k, local_dependence = spec$pairs,
                      variance_model = spec$variance, n_starts = n_starts, seed = seed,
                      tol = 1e-10, max_iter = 5000,
                      missing = if (isTRUE(spec$missing)) "fiml" else "error")
    if (!is.null(spec$min_variance)) arguments$min_variance <- spec$min_variance
    withCallingHandlers(do.call(lpa, arguments),
                        warning = function(w) {
                          # Bound-active fits are expected only where the case
                          # sets a floor; elsewhere the warning is kept.
                          if (!is.null(spec$floor) && inherits(w, "latents_boundary")) invokeRestart("muffleWarning")
                        })
  })
  scores <- vapply(fits, \(f) f$log_likelihood, numeric(1))
  best <- fits[[which.max(scores)]]
  best$seed_scores <- scores
  best
}

#' latents' class-level estimates and standard errors in the tidy layout
pc_latents_table <- function(fit, inference = TRUE) {
  k <- fit$n_profiles
  vars <- colnames(fit$means)
  estimates <- if (inference) parameter_inference(fit) else NULL
  se_of <- function(outcome, term, parameter) {
    if (is.null(estimates)) return(NA_real_)
    hit <- estimates$standard_error[estimates$outcome == outcome &
                                      estimates$term %in% term &
                                      estimates$parameter == parameter]
    if (length(hit) == 0L) NA_real_ else hit[1L]
  }
  pairs <- .latents_pattern_pairs(fit$covariance_pattern)
  equal <- identical(fit$variance_model, "equal")
  rows <- lapply(seq_len(k), function(j) {
    profile <- sprintf("profile_%d", j)
    owner <- if (equal) "shared" else profile
    sigma <- matrix(fit$covariances[, , j], length(vars), dimnames = list(vars, vars))
    split <- strsplit(pairs, ":", fixed = TRUE)
    rbind(
      data.frame(class = j, quantity = "proportion", term = "",
                 estimate = fit$profile_probabilities[1L, j],
                 standard_error = se_of(profile, "group_class_1", "probability")),
      data.frame(class = j, quantity = "mean", term = vars, estimate = fit$means[j, ],
                 standard_error = vapply(vars, \(v) se_of(profile, v, "mean"), numeric(1))),
      data.frame(class = j, quantity = "variance", term = vars, estimate = diag(sigma),
                 standard_error = vapply(vars, \(v) se_of(owner, sprintf("%s:%s", v, v),
                                                          "covariance"), numeric(1))),
      data.frame(class = j, quantity = "covariance", term = pairs,
                 estimate = vapply(split, \(p) sigma[p[1L], p[2L]], numeric(1)),
                 standard_error = vapply(split, \(p) se_of(owner, c(
                   sprintf("%s:%s", p[2L], p[1L]), sprintf("%s:%s", p[1L], p[2L])),
                   "covariance"), numeric(1))))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# -- Latent GOLD -------------------------------------------------------------------

pc_read_latin1 <- function(path) {
  sub("\r$", "", iconv(readLines(path, warn = FALSE), "latin1", "UTF-8"))
}

#' The last number on the first line matching `pattern`
pc_printed <- function(lines, pattern) {
  line <- grep(pattern, lines, value = TRUE)[1L]
  stopifnot("no line matches the pattern" = !is.na(line))
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  as.numeric(fields[length(fields)])
}

#' One row of the compact `write=` file, as numbers
pc_write_row <- function(lines, name) {
  line <- grep(sprintf(",%s,", name), lines, value = TRUE, fixed = TRUE)[1L]
  stopifnot("the write file has no such row" = !is.na(line))
  fields <- strsplit(line, ",", fixed = TRUE)[[1L]]
  as.numeric(fields[-seq_len(match(name, fields))])
}

#' Rows of a listing table (Variances, Covariances / Associations): the term
#' fields, the class, the printed coefficient and standard error.
pc_listing_rows <- function(lines, header, equal) {
  start <- which(lines == paste0(header, "\t"))[1L]
  stopifnot("the listing has no such table" = !is.na(start))
  body <- lines[(start + 2L):length(lines)]
  body <- body[seq_len(match(TRUE, !grepl("^y[0-9]", body)) - 1L)]
  fields <- strsplit(body, "\t", fixed = TRUE)
  do.call(rbind, lapply(fields, function(f) {
    cluster <- grep("^Cluster\\(", f, value = TRUE)
    numbers <- suppressWarnings(as.numeric(f))
    first_number <- which(!is.na(numbers))[1L]
    data.frame(first = f[1L], second = if (identical(f[2L], "<->")) f[3L] else f[1L],
               class = if (length(cluster)) as.integer(sub("Cluster\\((\\d+)\\)", "\\1", cluster))
                       else NA_integer_,
               printed = numbers[first_number], printed_se = numbers[first_number + 1L],
               printed_z = numbers[first_number + 2L])
  }))
}

#' Everything compared from one Latent GOLD run
pc_latentgold <- function(name, case) {
  spec <- case$spec
  k <- spec$k
  d <- length(spec$vars)
  equal <- identical(spec$variance, "equal")
  lines <- pc_read_latin1(pc_root("returned", paste0(name, ".lst")))
  written <- pc_read_latin1(pc_root("returned", paste0(name, "_write.txt")))
  parameters <- pc_write_row(written, "Parameters")
  parameter_se <- pc_write_row(written, "Parameters-SE")
  profile <- pc_write_row(written, "Profile")
  profile_se <- pc_write_row(written, "Profile-SE")
  npar <- as.integer(pc_printed(lines, "^Number of parameters \\(Npar\\)"))
  stopifnot("Parameters row length differs from Npar" = length(parameters) == npar)
  # Profile: k sizes, then per indicator the k class means and the overall one
  # (the overall means follow all class means when there are several).
  sizes <- profile[seq_len(k)]
  means <- matrix(profile[k + seq_len(k * d)], k, d)
  mean_se <- matrix(profile_se[k + seq_len(k * d)], k, d)
  # Variances and covariances: the full-precision Parameters entry whose value,
  # standard error and z-value all round to what the listing prints (the
  # z-value tells apart variances that print as 0.0000). The write file
  # orders them by block, not by the listing, so they are matched, not read
  # by position; a match must be unique.
  match_entry <- function(printed, printed_se, printed_z) {
    hit <- which(abs(parameters - printed) <= 5.0001e-5 &
                   abs(parameter_se - printed_se) <= 5.0001e-5 &
                   abs(parameters / parameter_se - printed_z) <= 5.0001e-5)
    stopifnot("a printed parameter matches no unique written one" = length(hit) == 1L)
    hit
  }
  covariance_rows <- function(header) {
    rows <- pc_listing_rows(lines, header, equal)
    rows$index <- mapply(match_entry, rows$printed, rows$printed_se, rows$printed_z)
    rows
  }
  variances <- covariance_rows("Variances")
  covariances <- covariance_rows("Covariances / Associations")
  expand <- function(rows, quantity) {
    do.call(rbind, lapply(seq_len(k), function(j) {
      own <- if (equal) rows else rows[rows$class == j, , drop = FALSE]
      data.frame(class = j, quantity = quantity,
                 term = if (quantity == "variance") own$first else
                   paste(own$first, own$second, sep = ":"),
                 estimate = parameters[own$index], standard_error = parameter_se[own$index])
    }))
  }
  table <- rbind(
    data.frame(class = seq_len(k), quantity = "proportion", term = "", estimate = sizes,
               standard_error = profile_se[seq_len(k)]),
    data.frame(class = rep(seq_len(k), d), quantity = "mean",
               term = rep(spec$vars, each = k), estimate = as.vector(means),
               standard_error = as.vector(mean_se)),
    expand(variances, "variance"), expand(covariances, "covariance"))
  saved <- utils::read.delim(pc_root("returned", paste0(name, "_posteriors.txt")),
                             check.names = FALSE)
  list(log_likelihood = pc_write_row(written, "Statistics")[1L],
       log_likelihood_printed = pc_printed(lines, "^Log-likelihood \\(LL\\)"),
       n_parameters = npar, bic = pc_printed(lines, "^BIC \\(based on LL\\)"),
       table = table, id = as.integer(saved$id),
       posteriors = as.matrix(saved[paste0("Cluster#", seq_len(k))]))
}

# -- Mplus -------------------------------------------------------------------------

pc_read_lines <- function(path) readLines(path, warn = FALSE)

#' TECH1 parameter numbers per class: NU (means) and THETA (lower triangle)
pc_mplus_tech1 <- function(lines, k, vars) {
  d <- length(vars)
  lapply(seq_len(k), function(j) {
    start <- grep(sprintf("PARAMETER SPECIFICATION FOR LATENT CLASS %d$", j), lines)
    stopifnot("no TECH1 block for a class" = length(start) == 1L)
    block <- lines[start:(start + 40L)]
    nu_at <- grep("^\\s+NU\\s*$", block)[1L]
    nu <- as.integer(strsplit(trimws(block[nu_at + 3L]), "\\s+")[[1L]])
    theta_at <- grep("^\\s+THETA\\s*$", block)[1L]
    theta <- matrix(0L, d, d, dimnames = list(vars, vars))
    invisible(lapply(seq_len(d), function(i) {
      values <- strsplit(trimws(block[theta_at + 2L + i]), "\\s+")[[1L]]
      stopifnot("THETA row label differs" = identical(tolower(values[1L]), vars[i]))
      theta[i, seq_len(i)] <<- as.integer(values[-1L])
    }))
    list(nu = nu, theta = theta)
  })
}

#' Mplus's ALPHA(C) parameter numbers (class logits, last class fixed at 0)
pc_mplus_alpha <- function(lines, k) {
  at <- grep("^\\s+ALPHA\\(C\\)\\s*$", lines)[1L]
  as.integer(strsplit(trimws(lines[at + 3L]), "\\s+")[[1L]])
}

#' Column order of a Mplus save file, from the output's SAVEDATA section
pc_mplus_save_columns <- function(lines) {
  start <- grep("Order and format of variables", lines)[1L]
  stop <- grep("Save file format", lines)
  stop <- stop[stop > start][1L]
  fields <- trimws(lines[(start + 1L):(stop - 1L)])
  fields <- fields[nzchar(fields)]
  toupper(vapply(strsplit(fields, "\\s+"), `[[`, character(1), 1L))
}

#' Everything compared from one Mplus run
pc_mplus <- function(name, case) {
  spec <- case$spec
  k <- spec$k
  vars <- spec$vars
  lines <- pc_read_lines(pc_root("returned", paste0(name, ".out")))
  npar <- as.integer(pc_printed(lines, "^Number of Free Parameters"))
  values <- scan(pc_root("returned", paste0(name, "-results.dat")), quiet = TRUE)
  tech1 <- pc_mplus_tech1(lines, k, vars)
  alpha <- pc_mplus_alpha(lines, k)
  # The RESULTS file holds every TECH1 parameter (an equality in MODEL
  # CONSTRAINT makes one of them dependent, so there can be more than the
  # free ones), then their standard errors, then the fit block.
  written <- max(unlist(lapply(tech1, \(t) c(t$nu, t$theta))), alpha)
  stopifnot("results file does not end with the fit block" =
              isTRUE(all.equal(values[2L * written + 1L], npar)))
  estimate <- \(i) ifelse(i > 0L, values[pmax(i, 1L)], NA_real_)
  se <- \(i) ifelse(i > 0L, values[written + pmax(i, 1L)], NA_real_)
  npar_offset <- written
  logits <- c(estimate(alpha[seq_len(k - 1L)]), 0)
  proportions <- exp(logits - max(logits)) / sum(exp(logits - max(logits)))
  split <- strsplit(spec$pairs, ":", fixed = TRUE)
  table <- do.call(rbind, lapply(seq_len(k), function(j) {
    theta <- tech1[[j]]$theta
    pair_index <- vapply(split, \(p) max(theta[p[1L], p[2L]], theta[p[2L], p[1L]]), integer(1))
    rbind(
      data.frame(class = j, quantity = "proportion", term = "", estimate = proportions[j],
                 standard_error = NA_real_),
      data.frame(class = j, quantity = "mean", term = vars, estimate = estimate(tech1[[j]]$nu),
                 standard_error = se(tech1[[j]]$nu)),
      data.frame(class = j, quantity = "variance", term = vars,
                 estimate = estimate(diag(theta)), standard_error = se(diag(theta))),
      data.frame(class = j, quantity = "covariance", term = spec$pairs,
                 estimate = estimate(pair_index), standard_error = se(pair_index)))
  }))
  saved <- utils::read.table(pc_root("returned", paste0(name, "-posteriors.dat")),
                             na.strings = "*")
  columns <- pc_mplus_save_columns(lines)
  stopifnot("save file does not match its layout" = length(columns) == ncol(saved))
  names(saved) <- columns
  list(log_likelihood = values[2L * npar_offset + 2L], n_parameters = npar,
       bic = values[2L * npar_offset + 4L], table = table, id = as.integer(saved$ID),
       posteriors = as.matrix(saved[grep("^CPROB", columns)]),
       terminated = any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", lines)),
       replicated = any(grepl("THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED", lines)),
       warnings = grep("WARNING", lines, value = TRUE))
}

# -- comparison -----------------------------------------------------------------------

#' latents' table relabelled into the reference's classes and joined with it
pc_join <- function(mine, theirs, order) {
  # order[j] is latents' class that matches the reference's class j.
  relabel <- match(mine$class, order)
  mine$class <- relabel
  merged <- merge(theirs, mine, by = c("class", "quantity", "term"),
                  suffixes = c("_reference", "_latents"))
  stopifnot("every reference quantity has a latents counterpart" =
              nrow(merged) == nrow(theirs))
  merged[order(merged$quantity, merged$class, merged$term), ]
}

#' Posterior matrix of latents in the reference's row order
pc_posteriors <- function(fit, id) {
  posteriors <- fit$subject_posteriors
  posteriors[id, , drop = FALSE]
}
