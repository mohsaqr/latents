suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
options(digits = 5, width = 160)
activities <- c("time_with_friends", "on_social_media", "tv_video_games","listened_music", "sports", "walking", "reading","part_time_job")
## A. workflow-lca filter: "p_adjusted" of get_results vs of diagnostics
set.seed(1)
f3 <- lca(student_esm, activities, n_classes = 3, n_starts = 5, seed = 1)
a <- get_results(f3, "residuals"); b <- suppressWarnings(get_residuals(f3))
cat("tests:", sum(!is.na(a$p_value)), " get_results p_adjusted<.01:", sum(a$p_adjusted < .01, na.rm=TRUE), " get_residuals p_adjusted<.01:", sum(b$p_adjusted < .01, na.rm=TRUE), " raw==adj:", isTRUE(all.equal(a$p_value,a$p_adjusted)), "\n")
print(subset(a, p_adjusted < .01)[c("profile","indicator_1","indicator_2","statistic","p_value","p_adjusted")])
print(subset(b, p_adjusted < .01)[c("profile","indicator_1","indicator_2","statistic","p_value","p_adjusted")])
print(attributes(a)[setdiff(names(attributes(a)), c("names","row.names","class"))])
## B. BLRT helper with invalid replicates: p-value not withheld
rep <- data.frame(valid = c(rep(TRUE, 15), rep(FALSE, 4)), statistic = c(rep(1, 14), 50, rep(NA, 4)))
print(suppressWarnings(latents:::.latents_blrt_p_value(rep, observed = 10, iter = 19)))
## C. printed table: big intercept + tiny slope
set.seed(2)
n <- 600; x <- runif(n, 0, 50000); cl <- rbinom(n, 1, .5)
d <- data.frame(y = ifelse(cl == 1, 200 + 0.004 * x, 500 - 0.003 * x) + rnorm(n, 0, 5), x = x)
mr <- mixture_regression(y ~ x, d, n_classes = 2, n_starts = 3, seed = 1)
tab <- get_results(mr, "coefficients"); print(tab); print(as.data.frame(tab)[c("class","term","estimate","std_error","conf_low","conf_high")])
## D. truth with n_group_classes = 1 and a group-level truth
f1 <- multilpa(course_engagement, c("browse","lectures","forum_read"), id = "student", n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
print(get_results(f1, "assignments", data = course_engagement, truth = "student_type"))
## E. plot default for transition fit
vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
mv <- lta(course_engagement, vars, "student", n_profiles = 3, n_group_classes = 2, time = "sequence", n_starts = 2, seed = 1)
p <- plot(mv); print(p$labels$title)
## prevalence by hand: expected share of class's enrolments in profile (posterior)
it <- get_results(mv, "initial")
print(it)
pp <- get_results(mv, "posteriors", format = "wide"); 
gp <- get_results(mv, "group_posteriors", format = "wide")
print(colMeans(as.matrix(pp[grep("^posterior_profile", names(pp))])))
A <- mv$transition_probabilities; delta <- mv$initial_probabilities %||% NULL; print(names(mv)[grepl("prob|prevalence", names(mv))]);
## F. occasion labels T1 ... T10
d10 <- transform(subset(course_engagement, sequence <= 10), wave = paste0("T", sequence))
print(tryCatch({lta(d10, vars[1:2], "student", n_profiles = 2, time = "wave", n_starts = 1, seed = 1); "fitted"}, error = function(e) paste("ERR", class(e)[1])))
## G. grid summary claims
grid <- enumerate_classes(course_engagement, vars[1:3], "student", n_profiles = 2:4, n_group_classes = 1:2, model = "VVI", n_starts = 2, seed = 1)
s <- summary(grid); print(s)
## H. tree masses by hand
tr <- get_tree(grid); print(head(tr, 4))
a2 <- get_candidate(grid, n_profiles = 2, n_group_classes = 1); a3 <- get_candidate(grid, n_profiles = 3, n_group_classes = 1)
m <- crossprod(a2$subject_posteriors, a3$subject_posteriors) / nrow(a2$subject_posteriors)
one <- subset(tr, from_candidate == "profiles2_groups1" & to_candidate == "profiles3_groups1")
print(one); print(m); print(rowSums(m)); print(colSums(m))
## I. ICC by hand, unbalanced, via anova; and within a stratum
set.seed(5); g <- rep(1:30, times = sample(2:12, 30, TRUE)); y <- rnorm(30)[g] + rnorm(length(g))
aov1 <- anova(lm(y ~ factor(g))); msb <- aov1[1,3]; msw <- aov1[2,3]; n <- length(y); k <- 30; ni <- tabulate(g); n0 <- (n - sum(ni^2)/n)/(k-1)
print(c(hand = (msb-msw)/(msb+(n0-1)*msw), pkg = descriptives(data.frame(y=y,g=g), vars="y", id="g")$icc))
dd <- data.frame(y=y,g=g,s=rep(c("a","b"),length.out=length(y)))
print(descriptives(dd, vars="y", id="g", by="s"))
## J. descriptives on a weighted fit with zero weights
w <- transform(course_engagement, wt = ifelse(student <= 5, 0, 1))
wf <- multilpa(w, vars[1:3], id = "student", n_profiles = 2, n_group_classes = 2, weights = "wt", n_starts = 2, seed = 1)
print(descriptives(wf)); print(descriptives(wf, by = "profile")); print(nobs(wf)); print(wf$n_groups); print(subset(get_results(wf,"counts"), level=="groups"))
print(get_classification(wf)[c("level","class","n_modal")])
