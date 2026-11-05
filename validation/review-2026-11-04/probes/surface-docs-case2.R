suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
options(digits = 5, width = 140)
vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
candidates <- enumerate_classes(course_engagement, vars, id = "student",
                                n_profiles = 2:3, n_group_classes = 2:3,
                                model = c("VVI", "VVV"), n_starts = 6,
                                max_iter = 2000, tol = 1e-10, seed = 1)
tab <- get_results(candidates, "candidates")
print(tab[c("candidate","n_profiles","n_group_classes","model","log_likelihood","bic_groups","converged","boundary","n_best_replicated")])
full <- get_candidate(candidates, n_profiles = 2, n_group_classes = 3, model = "VVV")
print(head(get_results(full, "residuals"), 4))
covariate_fit <- multilpa(course_engagement, vars, id = "student", n_profiles = 3, n_group_classes = 3,
  profile_covariates = "previous_grade", n_starts = 6, max_iter = 2000, tol = 1e-10, seed = 1)
print(covariate_fit)
print(get_results(covariate_fit, "coefficients"))
staged <- fit_staged(course_engagement, vars, id = "student", n_profiles = 3, n_group_classes = 3,
  n_starts = 6, max_iter = 2000, tol = 1e-10, seed = 1)
print(staged)
print(get_results(staged, "stages"))
moves <- lta(course_engagement, vars, "student", n_profiles = 3, n_group_classes = 2, time = "sequence", n_starts = 10, seed = 1)
print(get_results(moves, "starts"))
print(get_results(moves))
print(get_results(moves, "initial"))
print(get_results(moves, "transitions"))
print(get_results(moves, "entropy"))
print(get_results(moves, "sequence_summary"))
r <- tryCatch(parameter_inference(moves), error = function(e) {cat("PI err:", class(e), conditionMessage(e), "\n"); NULL})
