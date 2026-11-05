suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
options(digits = 5, width = 170)
vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
fit <- multilpa(course_engagement, vars = vars, id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1)
direct <- multilpa(course_engagement, vars = vars, id = "student", n_profiles = 2, n_group_classes = 2,
                   direct_effects = c("attendance ~ previous_grade", "forum_post ~ previous_grade"), n_starts = 3, seed = 1)
print(get_results(direct, "direct_effects")); print(compare_models(without = fit, direct = direct))
print(get_results(fit)[c("profile","indicator","mean")]); print(get_results(direct)[c("profile","indicator","mean")])
# tree-plot two-class solution
activities <- c("time_with_friends", "on_social_media", "tv_video_games","listened_music", "sports", "walking", "reading","part_time_job")
m <- multilpa(student_esm, activities, "student", n_profiles = 2, n_group_classes = 1, categorical = activities, n_starts = 3, seed = 1)
print(reshape(subset(get_results(m,"responses"), category=="yes")[c("profile","indicator","probability")], idvar="indicator", timevar="profile", direction="wide")); print(subset(get_results(m,"counts"), level=="individuals"))
print(sapply(student_esm[activities], function(v) mean(v == "yes")))
# case-esm: "Profile 1 is therefore media-heavy"; class numbering in the 2x2
# evaluation reference
s <- sensitivity(fit, seeds = 1:2); print(s)
# identical on deleted rows at a tight tolerance
set.seed(3)
empty <- sample(nrow(course_engagement), 150); partial <- setdiff(sample(nrow(course_engagement), 150), empty)
inc <- transform(course_engagement, browse = replace(browse, empty, NA), lectures = replace(lectures, empty, NA), forum_read = replace(forum_read, empty, NA), forum_post = replace(forum_post, empty, NA), attendance = replace(attendance, c(empty, partial), NA))
a <- multilpa(inc, vars, "student", n_profiles = 3, n_group_classes = 2, model = "VVI", missing = "fiml", n_starts = 3, tol = 1e-10, max_iter = 3000, seed = 1)
d <- multilpa(inc[-empty, ], vars, "student", n_profiles = 3, n_group_classes = 2, model = "VVI", missing = "fiml", n_starts = 3, tol = 1e-10, max_iter = 3000, seed = 1)
print(c(a$log_likelihood, d$log_likelihood, a$bic_individual, d$bic_individual)); print(get_entropy(a)$relative_entropy); print(get_entropy(d)$relative_entropy)
