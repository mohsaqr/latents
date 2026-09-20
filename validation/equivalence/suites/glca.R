# Nonparametric two-level latent class analysis against glca (Kim, Jeon, Chang
# & Chung 2022, Applied Psychological Measurement 46(5):439-441).
#
# glca is an independent implementation of the same model: its documentation
# states that the "level-2 (group-level) latent class is categorized by the
# prevalence of level-1 latent class for group variable", and it cites Vermunt
# (2003) directly. Its `measure.inv = TRUE` default is multilpa's shared
# measurement model. The parameter blocks map one to one:
#
#   glca delta  <-> multilpa group_probabilities     P(group class)
#   glca gamma  <-> multilpa profile_probabilities   P(profile | group class)
#   glca rho    <-> multilpa response_probabilities  P(item response | profile)
#
# Both datasets carry missing item responses, so this also checks that the two
# observed-data likelihoods agree, not only the complete-data ones.

suite_glca <- function() {
  require_suite_packages(
    "glca",
    reason = "the independent Vermunt (2003) two-level implementation, and the nyts18 and gss08 data, that this suite compares against")
  cases <- list(
    list(dataset = "nyts18", id = "SCH_ID",
         items = c("ECIGT", "ECIGAR", "ESLT", "EELCIGT", "EHOOKAH"),
         profiles = 3L, group_classes = 2L),
    list(dataset = "nyts18", id = "SCH_ID",
         items = c("ECIGT", "ECIGAR", "ESLT", "EELCIGT", "EHOOKAH"),
         profiles = 2L, group_classes = 2L),
    list(dataset = "gss08", id = "REGION",
         items = c("DEFECT", "HLTH", "RAPE", "POOR", "SINGLE", "NOMORE"),
         profiles = 3L, group_classes = 2L),
    list(dataset = "gss08", id = "REGION",
         items = c("DEFECT", "HLTH", "RAPE", "POOR", "SINGLE", "NOMORE"),
         profiles = 2L, group_classes = 3L)
  )
  do.call(rbind, lapply(cases, function(case) {
    do.call(.glca_case, case)
  }))
}

#' Run one glca comparison
#' @param dataset Name of a glca dataset.
#' @param id Name of the level-2 identifier column.
#' @param items Character vector of binary item names.
#' @param profiles Number of level-1 latent classes.
#' @param group_classes Number of level-2 latent classes.
#' @return A `data.frame` of compared quantities.
.glca_case <- function(dataset, id, items, profiles, group_classes) {
  observed <- get(utils::data(list = dataset, package = "glca",
                              envir = environment()), envir = environment())
  formula <- stats::as.formula(sprintf("glca::item(%s) ~ 1", paste(items, collapse = ", ")))
  reference <- glca::glca(formula, group = observed[[id]], data = observed,
                          nclass = profiles, ncluster = group_classes,
                          n.init = 20L, seed = 1L, maxiter = 20000L,
                          eps = 1e-10, verbose = FALSE)
  fit <- multilpa(observed, vars = items, id = id,
                  n_profiles = profiles, n_group_classes = group_classes,
                  categorical = items, n_starts = 30L, seed = 1L,
                  tol = 1e-12, max_iter = 20000L, missing = "fiml")

  profile_order <- align_classes(flatten_blocks(reference$param$rho),
                                 flatten_blocks(fit$response_probabilities))
  responses <- lapply(fit$response_probabilities,
                      \(block) unname(as.matrix(block))[profile_order, , drop = FALSE])
  prevalence <- unname(fit$profile_probabilities)[, profile_order, drop = FALSE]
  group_order <- align_classes(unname(reference$param$gamma), prevalence)
  prevalence <- prevalence[group_order, , drop = FALSE]

  label <- sprintf("%s, %d profiles x %d group classes", dataset, profiles, group_classes)
  n_categories <- ncol(reference$param$rho[[1L]])
  cell_names <- sprintf("%s: P(%s = level %d | profile %d)", label,
                        rep(items, each = profiles * n_categories),
                        rep(rep(seq_len(n_categories), each = profiles), length(items)),
                        rep(seq_len(profiles), n_categories * length(items)))

  rbind(
    compare_values(sprintf("%s: maximised log-likelihood", label),
                   reference$gof$loglik, as.numeric(logLik(fit)), tolerance = 1e-5),
    compare_values(sprintf("%s: AIC", label), reference$gof$AIC, fit$aic, tolerance = 1e-4),
    compare_values(sprintf("%s: BIC", label), reference$gof$BIC, fit$bic_individual,
                   tolerance = 1e-4),
    compare_values(sprintf("%s: group-class probability %d", label, seq_len(group_classes)),
                   sort(as.numeric(reference$param$delta)),
                   sort(as.numeric(fit$group_probabilities)), tolerance = 1e-4),
    compare_values(sprintf("%s: P(profile %d | group class %d)", label,
                           rep(seq_len(profiles), each = group_classes),
                           rep(seq_len(group_classes), profiles)),
                   as.numeric(unname(reference$param$gamma)),
                   as.numeric(prevalence), tolerance = 1e-4),
    compare_values(cell_names,
                   as.numeric(vapply(reference$param$rho, \(b) unname(as.matrix(b)),
                                     matrix(0, profiles, n_categories))),
                   as.numeric(vapply(responses, identity, matrix(0, profiles, n_categories))),
                   tolerance = 1e-4)
  ) |> transform(source = sprintf("glca %s", utils::packageVersion("glca")),
                 dataset = dataset)
}
