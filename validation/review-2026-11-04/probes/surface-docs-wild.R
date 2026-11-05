suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
options(digits = 5, width = 170)
vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
check <- function(fit, name) {
  cat("\n=====", name, "\n")
  tab <- tryCatch(get_classification(fit), error = function(e) {cat("get_classification ERR", class(e)[1], conditionMessage(e), "\n"); NULL})
  ent <- tryCatch(get_entropy(fit), error = function(e) {cat("get_entropy ERR", class(e)[1], "\n"); NULL})
  ave <- tryCatch(get_avepp(fit), error = function(e) {cat("get_avepp ERR", class(e)[1], "\n"); NULL})
  if (!is.null(tab)) print(tab[c("level","class","n_modal","estimated_n","estimated_proportion","average_posterior")])
  if (!is.null(ent)) print(ent)
  pick <- function(p, lvl) if (inherits(p, "latents_plots")) p[[lvl]] else p
  for (lvl in unique(tab$level)) {
    t <- tab[tab$level == lvl, ]
    s <- pick(plot_sizes(fit), lvl)$data
    cat(lvl, "sizes equal:", isTRUE(all.equal(unname(s$count), t$estimated_n, tolerance = 1e-8)), " shares equal:", isTRUE(all.equal(unname(s$share), t$estimated_proportion, tolerance=1e-8)), "\n")
    pp <- pick(plot_posteriors(fit), lvl)
    d <- pp$data; cl <- if ("class" %in% names(d)) d$class else d$profile
    cat(lvl, "n points:", nrow(d), " modal counts plot:", paste(tabulate(cl), collapse=","), " table:", paste(t$n_modal, collapse=","), "\n")
    pe <- pick(plot_entropy(fit), lvl); cat(lvl, "subtitle:", pe$labels$subtitle, "\n")
    a <- pick(plot_avepp(fit), lvl)$data; a <- a[a$row == a$column, ]; a <- a[order(a$row),]
    cat(lvl, "avepp diag plot:", paste(round(a$average_posterior, 5), collapse=","), " table:", paste(round(t$average_posterior, 5), collapse=","), "\n")
  }
}
gr2 <- mixture_regression(score ~ wave, growth_schools, n_classes = 2, id = "student", class_level = "group", random = "wave", random_covariance = "equal", cluster = "school", n_group_classes = 2, n_starts = 2, seed = 1)
check(gr2, "multilevel growth")
reg2 <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2, id = "student", n_group_classes = 2, n_starts = 3, seed = 1); check(reg2, "two-level regression")
print(get_results(reg2, "classes")); print(get_results(reg2, "group_classes"))
add <- multilpa(course_engagement, vars[1:3], "student", n_group_classes = 2, family = "additive", seed = 1); check(add, "additive")
fc <- multilpa(course_engagement, vars[1:3], "student", n_profiles = 2, n_group_classes = 2, family = "full_cross_level", seed = 1); check(fc, "full cross level")
print(get_results(fc,"composition"))
gl <- lta(course_engagement, vars[1:3], "student", n_profiles = 2, n_group_classes = 2, time = "sequence", transitions = "occasion", n_starts = 2, seed = 1, max_iter = 300); check(gl, "general lta")
wr <- transform(study_hours, wt = ifelse(student %% 3 == 0, 2, 1))
wreg <- mixture_regression(score ~ hours, data = wr, n_classes = 2, id = "student", n_group_classes = 2, weights = "wt", n_starts = 3, seed = 1); check(wreg, "weighted two-level regression")
## grid views vs tables
grid <- enumerate_classes(course_engagement, vars[1:3], "student", n_profiles = 1:3, n_group_classes = 1:2, model = "VVI", n_starts = 2, seed = 1)
pc <- plot_criteria(grid); print(pc$labels$subtitle)
pq <- plot_entropy(grid); print(pq$labels$subtitle)
dl <- plot_deltas(grid); print(head(dl$data, 8)); print(get_deltas(grid)[c("candidate","delta_profile","delta_group_class")])
print(subset(as.data.frame(grid), TRUE)[c("candidate","converged","error")])
sm <- plot_sizes(grid); print(sm$labels$subtitle); print(sm$data[c("candidate","smallest_share")])
ct <- get_classification(grid); print(head(ct, 3))
hand <- with(subset(ct, level == "individuals"), tapply(estimated_proportion, candidate, min)); print(hand)
tr <- plot_trust(grid); print(tr$data[c("candidate","status")])
# staged vs joint, test = lrt
f0 <- multilpa(course_engagement, vars, id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1)
st <- fit_staged(course_engagement, vars, id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1)
print(tryCatch(compare_models(joint = f0, staged = st, test = "lrt"), error = function(e) paste("ERR", class(e)[1], conditionMessage(e))))
cmp <- compare_models(joint = f0, staged = st); print(as.data.frame(cmp)[c("model","n_parameters","log_likelihood","aic","bic_groups","delta_bic","bic_weight")])
print(c(hand_with_measurement = -2 * st$log_likelihood + st$n_parameters_with_measurement * log(106), joint = f0$bic))
print(get_results(st, "information_criteria"))
big <- enumerate_classes(course_engagement, vars, id = "student", n_profiles = 3, n_group_classes = 2, model = "VVI", n_starts = 2, seed = 1)
print(tryCatch(class(bootstrap_lrt(f0, get_candidate(big, n_profiles = 3, n_group_classes = 2), iter = 2, seed = 42)), error = function(e) paste("ERR", class(e)[1], conditionMessage(e))))
