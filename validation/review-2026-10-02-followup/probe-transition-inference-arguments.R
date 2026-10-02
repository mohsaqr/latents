devtools::load_all(".", quiet = TRUE)
vars <- c("browse", "lectures")
vars <- intersect(vars, names(course_engagement))
if (length(vars) < 2) vars <- names(course_engagement)[vapply(course_engagement, is.numeric, logical(1))][3:4]
cat("vars:", vars, "\n")
probe <- function(fit, label) {
  base <- parameter_inference(fit)
  for (args in list(list(method = "bootstrap", iter = 3, seed = 1), list(adjust = "bonferroni"),
                    list(boundary = "fix"))) {
    r <- tryCatch(do.call(parameter_inference, c(list(fit), args)), error = \(e) e)
    cat(sprintf("%-12s %-30s -> %s\n", label, paste(names(args), args, sep = "=", collapse = ","),
      if (inherits(r, "error")) paste("ERROR", class(r)[1]) else
        if (identical(r, base)) "IDENTICAL to Wald default (silently ignored)" else "different result"))
  }
}
hom <- lta(course_engagement, vars, "student", n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
gen <- lta(course_engagement, vars, "student", n_profiles = 2, time = "sequence",
           transition_covariates = "previous_grade", n_starts = 2, seed = 1)
cat(class(hom)[1], class(gen)[1], "\n")
probe(hom, class(hom)[1]); probe(gen, class(gen)[1])
