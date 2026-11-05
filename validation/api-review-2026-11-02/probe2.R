suppressMessages(devtools::load_all(".", quiet=TRUE))
try_msg <- function(expr) tryCatch({force(expr); "OK"}, error=function(e) paste0("[", class(e)[1], "] ", conditionMessage(e)), warning=function(w) paste0("[W ", class(w)[1], "] ", conditionMessage(w)))
vars <- c("browse","lectures","forum_read","forum_post","attendance")
fits <- list(
  lpa = lpa(faithful, c("eruptions","waiting"), n_profiles=2, n_starts=2, seed=1),
  lca = multilca(subset(student_esm, day <= 1), vars=c("time_with_friends","on_social_media","tv_video_games"), id="student", n_profiles=2, n_group_classes=2, n_starts=1, seed=1),
  lta = lta(course_engagement, vars=vars, id="student", n_profiles=2, time="sequence", n_starts=2, seed=1),
  mixreg = mixture_regression(score ~ hours, data=study_hours, n_classes=2, n_starts=3, seed=1),
  growth = mixture_regression(score ~ wave, subset(growth_scores, student <= 60), n_classes=2, id="student", class_level="group", random="intercept", n_starts=3, seed=1)
)
grid <- enumerate_classes(subset(course_engagement, student <= 30), c("browse","lectures"), "student", n_profiles=2:3, n_group_classes=1:2, model="VVI", n_starts=2, seed=1)
for (n in names(fits)) {   # probe loop over 5 fits, diagnostic script only
  f <- fits[[n]]
  cat("\n=====", n, "class:", class(f)[1], "\n")
  cat(" report:", try_msg(capture.output(report(f, plots=FALSE))), "\n")
  cat(" diagnostics:", try_msg(capture.output(d <- diagnostics(f, plots=FALSE))), "\n")
  if (exists("d")) { cat(" diag tables:", try_msg(get_results(d, "zz")), "\n"); rm(d) }
  cat(" plot what=zz:", try_msg(plot(f, what="zz")), "\n")
  cat(" plot_diagnostics:", try_msg(plot_diagnostics(f)), "\n")
  cat(" plot_profiles:", try_msg(plot_profiles(f)), "\n")
  cat(" parameter_inference:", try_msg(capture.output(parameter_inference(f))), "\n")
  cat(" descriptives:", try_msg(descriptives(f)), "\n")
}
cat("\n===== grid\n")
cat(" report:", try_msg(capture.output(report(grid, plots=FALSE))), "\n")
cat(" diagnostics:", try_msg(capture.output(diagnostics(grid, plots=FALSE))), "\n")
cat(" plot what=zz:", try_msg(plot(grid, what="zz")), "\n")
cat(" plot_views(grid):\n"); print(plot_views(grid))
cat(" plot_enumeration criterion:", try_msg(plot_enumeration(grid, what="bic")), "\n")
cat(" plot(grid, criterion='bic'):", try_msg(plot(grid, criterion="bic")), "\n")
cat(" plot_selection:", try_msg(plot_selection(grid)), "\n")
