# Helpers of the ML three-step feature kit: fitting a case, writing the
# Latent GOLD and Mplus inputs, and reading their outputs. Sourced by
# make-kit.R, make-pipeline-step3.R and equivalence/test-features-ml-three-step.R
# (never a helper-*.R file, which testthat would load for every test file).

# From the project root (the kit scripts) or from equivalence/ (testthat).
ml3_root <- function() {
  here <- file.path("features", "ml-three-step")
  if (dir.exists(here)) here else file.path("equivalence", here)
}
ml3_kit <- function() file.path(ml3_root(), "kit")
ml3_returned <- function() file.path(ml3_root(), "returned")

# -- fitting --------------------------------------------------------------------

#' The latents measurement fit of a case (deterministic: fixed seed and starts)
ml3_fit <- function(case) {
  spec <- case$spec
  common <- list(n_starts = spec$n_starts %||% 20L, seed = 1L, tol = 1e-10, max_iter = 5000L)
  extra <- spec$extra %||% list()
  fit <- switch(spec$family,
    lpa = do.call(lpa, c(list(case$data, spec$vars, n_profiles = spec$k,
                              model = spec$model %||% "VVI"), common, extra)),
    lca = do.call(lca, c(list(case$data, spec$vars, n_classes = spec$k), common, extra)),
    ml_lpa = do.call(multilpa, c(list(case$data, spec$vars, "g", n_profiles = spec$k,
                                      n_group_classes = spec$h, model = spec$model %||% "VVI"),
                                 common, extra)))
  fit
}

#' Map fixed-seed cases without forking an unsafe BLAS runtime
#'
#' `ml3_fit()` sets its own seed for each case. Serial execution therefore
#' preserves the fitting stream and is also available explicitly through
#' `options(latents.ml_three_step.serial = TRUE)` for reproducibility.
#' @param cases A list of cases.
#' @param fun A function of one case.
#' @return A list in input order.
ml3_map_cases <- function(cases, fun) {
  if (isTRUE(getOption("latents.ml_three_step.serial", FALSE)) ||
      identical(.Platform$OS.type, "windows") || !latents:::.latents_fork_safe()) {
    lapply(cases, fun)
  } else {
    parallel::mclapply(cases, fun, mc.cores = 2L, mc.preschedule = FALSE)
  }
}

#' Posteriors, modal classes, clusters and outcome at the analysed level
ml3_level_inputs <- function(fit, case) {
  if (identical(case$level, "groups")) {
    groups <- split(case$data[[case$outcome]], fit$group_index)
    list(posteriors = unname(fit$group_posteriors), modal = fit$group_classes,
         values = vapply(groups, \(v) v[[1L]], numeric(1)),
         clusters = seq_along(fit$group_classes), ids = seq_along(fit$group_classes))
  } else {
    list(posteriors = unname(fit$subject_posteriors), modal = fit$subject_profiles,
         values = case$data[[case$outcome]], clusters = fit$group_index,
         ids = case$data$id)
  }
}

ml3_entropy <- function(posteriors) {
  if (ncol(posteriors) < 2L) return(NA_real_)
  1 - sum(-posteriors * log(pmax(posteriors, 1e-300))) /
    (nrow(posteriors) * log(ncol(posteriors)))
}

#' The permutation of `mine`'s columns closest to `theirs`
ml3_relabel <- function(mine, theirs) {
  permutations <- function(k) if (k == 1L) list(1L) else
    do.call(c, lapply(seq_len(k), \(i) lapply(permutations(k - 1L),
                                             \(p) c(i, setdiff(seq_len(k), i)[p]))))
  candidates <- permutations(ncol(mine))
  gaps <- vapply(candidates, \(p) max(abs(mine[, p, drop = FALSE] - theirs)), numeric(1))
  list(order = candidates[[which.min(gaps)]], gap = min(gaps))
}

# -- Latent GOLD -----------------------------------------------------------------

ml3_lg_models <- function(case) {
  grid <- expand.grid(se = c("robust", "standard"), variance = c("varying", "equal"),
                      stringsAsFactors = FALSE)
  if (isTRUE(case$psu)) {
    grid <- rbind(grid, data.frame(se = "psu", variance = c("varying", "equal")))
  }
  grid$tag <- paste0(substr(grid$variance, 1L, 1L),
                     c(robust = "r", standard = "s", psu = "p")[grid$se])
  grid
}

#' A step-three .lgs file reading the posterior columns of `data_file`
ml3_lg_step3 <- function(name, data_file, k, models, outcome = "y") {
  posterior <- paste(sprintf("Cluster#%d", seq_len(k)), collapse = " ")
  blocks <- lapply(seq_len(nrow(models)), \(i) {
    m <- models[i, ]
    c("model", sprintf("title '%s_%s';", name, m$tag), "options", "   maxthreads=1;",
      "   algorithm",
      "      tolerance=1e-012 emtolerance=1e-010 emiterations=10000 nriterations=1000;",
      "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
      "   missing excludeall;", "   step3 modal ml;", "   output",
      sprintf("      parameters=first standarderrors%s profile estimatedvalues=model",
              if (identical(m$se, "robust")) "=robust" else ""),
      sprintf("      write='%s_%s.txt';", name, m$tag),
      "variables",
      if (identical(m$se, "psu")) "   psuid g;",
      sprintf("   dependent %s continuous;", outcome),
      sprintf("   latent Cluster nominal posterior = ( %s );", posterior),
      "equations", sprintf("   %s <- 1 + Cluster;", outcome),
      if (identical(m$variance, "varying")) sprintf("   %s | Cluster;", outcome),
      "end model", "")
  })
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s'", data_file), "", unlist(blocks))
}

#' A step-three data file: id, g, outcome and full-precision posteriors
ml3_write_lg_data <- function(path, ids, clusters, values, posteriors, outcome = "y") {
  columns <- c(list(id = ids, g = clusters),
               stats::setNames(list(values), outcome),
               stats::setNames(lapply(seq_len(ncol(posteriors)), \(t) posteriors[, t]),
                               sprintf("Cluster#%d", seq_len(ncol(posteriors)))))
  text <- vapply(columns, \(v) if (is.double(v)) sprintf("%.17g", v) else as.character(v),
                 character(length(ids)))
  text <- matrix(text, nrow = length(ids))
  lines <- c(paste(names(columns), collapse = "\t"), apply(text, 1L, paste, collapse = "\t"))
  writeLines(lines, path, sep = "\r\n")
}

#' One Latent GOLD `write=` file: the fields of each output section
ml3_read_lg_write <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  fields <- strsplit(lines[-1L], ",", fixed = TRUE)
  named <- Filter(\(f) length(f) >= 5L && nzchar(f[5L]), fields)
  stats::setNames(lapply(named, \(f) as.numeric(f[-seq_len(5L)])),
                  vapply(named, `[`, character(1), 5L))
}

#' The step-three quantities of one Latent GOLD model
ml3_lg_result <- function(name, tag, k, variance) {
  path <- file.path(ml3_returned(), sprintf("%s_%s.txt", name, tag))
  if (!file.exists(path)) return(NULL)
  w <- ml3_read_lg_write(path)
  n_var <- if (identical(variance, "equal")) 1L else k
  pairs <- utils::combn(k, 2L)
  # WaldStatistics: Wald(0) intercept, Wald(0) Cluster, Wald(0) variances,
  # [Wald(=) variances], paired comparisons of the means, [of the variances].
  wald <- w$WaldStatistics
  head <- if (n_var > 1L) 4L else 3L
  list(log_likelihood = w$Statistics[1L],
       sizes = w$EstimatedValues[seq_len(k)],
       means = w$EstimatedValues[k + seq_len(k)],
       mean_se = w[["EstimatedValues-SE"]][k + seq_len(k)],
       variances = utils::tail(w$Parameters, n_var),
       variance_se = utils::tail(w[["Parameters-SE"]], n_var),
       contrasts_first = w$Parameters[1L + seq_len(k - 1L)],
       contrasts_first_se = w[["Parameters-SE"]][1L + seq_len(k - 1L)],
       pair_wald = wald[head + seq_len(ncol(pairs))],
       wald_head = wald[seq_len(head)])
}

#' A Latent GOLD step-one .lgs (pipeline): fit and write 6-decimal posteriors
ml3_lg_step1 <- function(name, data_file, spec) {
  categorical <- identical(spec$family, "lca")
  type <- if (categorical) "nominal" else "continuous"
  c("//LG6.1//", "version = 6.1", sprintf("infile '%s'", data_file), "",
    "model", sprintf("title '%s';", name), "options", "   maxthreads=1;",
    "   algorithm",
    "      tolerance=1e-012 emtolerance=1e-010 emiterations=10000 nriterations=1000;",
    "   startvalues", "      seed=20261020 sets=100 tolerance=1e-005 iterations=250;",
    "   bayes", "      categorical=0 variances=0 latent=0 poisson=0;",
    "   missing excludeall;", "   output", "      parameters=first profile;",
    sprintf("   outfile '%s_post.txt' classification keep id y;", name),
    "variables",
    sprintf("   dependent %s;", paste(sprintf("%s %s", spec$vars, type), collapse = ", ")),
    sprintf("   latent Cluster nominal %d;", spec$k),
    "equations", "   Cluster <- 1;", sprintf("   %s <- 1 + Cluster;", spec$vars),
    if (!categorical) sprintf("   %s | Cluster;", spec$vars), "end model")
}

ml3_read_lg_lines <- function(path) {
  sub("\r$", "", iconv(readLines(path, warn = FALSE), "latin1", "UTF-8"))
}

#' The step-one log likelihood of a Latent GOLD listing
ml3_lg_step1_loglik <- function(path) {
  lines <- ml3_read_lg_lines(path)
  line <- grep("^Log-likelihood \\(LL\\)", lines, value = TRUE)[1L]
  as.numeric(strsplit(trimws(line), "\t")[[1L]][2L])
}

# -- Mplus -----------------------------------------------------------------------

#' An Mplus step-one input with DU3STEP and DE3STEP for the same outcome
ml3_mplus_input <- function(name, spec, columns) {
  categorical <- identical(spec$family, "lca")
  vars <- paste(spec$vars, collapse = " ")
  class_block <- unlist(lapply(seq_len(spec$k), \(j) {
    if (categorical) return(c(sprintf(" %%c#%d%%", j),
                              sprintf(" [%s];", paste0(spec$vars, "$1", collapse = " "))))
    c(sprintf(" %%c#%d%%", j), sprintf(" [%s];", vars), sprintf(" %s;", vars))
  }))
  c(sprintf("TITLE: latents ML three-step %s;", name),
    sprintf("DATA: FILE = %s.dat;", name), "VARIABLE:",
    sprintf(" NAMES = %s;", paste(columns, collapse = " ")),
    sprintf(" USEVARIABLES = %s;", vars), " IDVARIABLE = id;",
    if (categorical) sprintf(" CATEGORICAL = %s;", vars),
    sprintf(" CLASSES = c(%d);", spec$k),
    " AUXILIARY = ydu (DU3STEP) yde (DE3STEP);",
    "ANALYSIS:", " TYPE = MIXTURE;", " STARTS = 200 40;", " STITERATIONS = 40;",
    " STSEED = 20261020;", " MITERATIONS = 10000;", " LOGCRITERION = 0.00000001;",
    " RLOGCRITERION = 0.000000001;", " MCONVERGENCE = 0.000001;",
    "MODEL:", " %OVERALL%", class_block,
    "SAVEDATA:", sprintf(" FILE = %s-post.dat;", name), " SAVE = CPROBABILITIES;",
    " FORMAT = F20.12;")
}

#' The step-one likelihood, saved posteriors and both auxiliary tables
ml3_mplus_result <- function(name, k) {
  path <- file.path(ml3_returned(), paste0(name, ".out"))
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  ll_line <- grep("^\\s*H0 Value\\s", lines, value = TRUE)[1L]
  saved <- as.matrix(utils::read.table(file.path(ml3_returned(), paste0(name, "-post.dat"))))
  # Column names from the output's SAVEDATA layout.
  first <- grep("Order and format of variables", lines)[1L]
  layout <- trimws(lines[first + seq_len(60L)])
  layout <- layout[seq_len(which(!nzchar(layout[-1L]))[1L])]
  layout <- layout[nzchar(layout)]
  colnames(saved) <- vapply(strsplit(layout, "\\s+"), `[`, character(1), 1L)
  posteriors <- saved[, sprintf("CPROB%d", seq_len(k)), drop = FALSE]
  start <- grep("^EQUALITY TESTS OF MEANS", lines)
  table_for <- function(variable) {
    at <- grep(sprintf("^%s\\s*$", toupper(variable)), lines)
    at <- at[at > start][1L]
    class_lines <- lines[at + seq_len(40L)]
    rows <- grep("^\\s*Class [0-9]+\\s+-?[0-9.]+\\s+[0-9.]+\\s*$", class_lines, value = TRUE)
    numbers <- t(vapply(strsplit(trimws(rows[seq_len(k)]), "\\s+"), \(f) as.numeric(f[3:4]),
                        numeric(2)))
    pair_rows <- grep("^\\s*Class [0-9]+ vs\\. [0-9]+", class_lines, value = TRUE)
    pair_fields <- strsplit(trimws(pair_rows), "\\s+")
    counts <- grep("Final Class Counts", class_lines)[1L]
    size_lines <- grep("^\\s*[0-9]+\\s+[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s*$",
                       class_lines[counts + seq_len(k + 8L)], value = TRUE)
    overall <- grep("^\\s*Overall test", class_lines, value = TRUE)[1L]
    list(means = numbers[, 1L], mean_se = numbers[, 2L],
         pairs = data.frame(
           class_a = vapply(pair_fields, \(f) as.integer(f[2L]), integer(1)),
           class_b = vapply(pair_fields, \(f) as.integer(f[4L]), integer(1)),
           chi_square = vapply(pair_fields, \(f) as.numeric(f[5L]), numeric(1))),
         overall = as.numeric(strsplit(trimws(overall), "\\s+")[[1L]][3L]),
         last_sizes = vapply(strsplit(trimws(size_lines), "\\s+"),
                             \(f) as.numeric(f[5L]), numeric(1)))
  }
  list(log_likelihood = as.numeric(strsplit(trimws(ll_line), "\\s+")[[1L]][3L]),
       posteriors = unname(posteriors), ids = saved[, "ID"],
       du3step = table_for("YDU"), de3step = table_for("YDE"))
}
