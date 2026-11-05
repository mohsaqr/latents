suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
options(digits = 6, width = 150)
library(ggplot2)
vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
ld <- function(p, i) ggplot2::layer_data(p, i)
## SE mapping for every covariance structure: profile-table SE == parameter_inference row
d <- subset(course_engagement, student <= 60)
models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","VVV","EEV","VEE","EVE","VVE","VEV","EVV")
res <- lapply(models, function(m) {
  f <- tryCatch(multilpa(d, vars[1:3], id = "student", n_profiles = 2, n_group_classes = 1, model = m, n_starts = 2, seed = 1), error = function(e) e)
  if (inherits(f, "error")) return(data.frame(model = m, status = paste("fit err", conditionMessage(f)), na_mean = NA, na_var = NA, mean_ok = NA))
  p <- suppressMessages(withCallingHandlers(get_results(f, "profiles"), warning = function(w) invokeRestart("muffleWarning")))
  pi <- tryCatch(as.data.frame(parameter_inference(f)), error = function(e) e)
  if (inherits(pi, "error")) return(data.frame(model = m, status = paste("PI err", class(pi)[1]), na_mean = sum(is.na(p$mean_standard_error)), na_var = sum(is.na(p$variance_standard_error)), mean_ok = NA))
  mp <- pi[pi$parameter == "mean", ]
  key <- paste0("profile_", p$profile, "|", p$indicator)
  se <- mp$standard_error[match(key, paste0(mp$outcome, "|", mp$term))]
  data.frame(model = m, status = paste(unique(pi$parameter[pi$level == "measurement"]), collapse = ","), na_mean = sum(is.na(p$mean_standard_error)), na_var = sum(is.na(p$variance_standard_error)), mean_ok = isTRUE(all.equal(se, p$mean_standard_error)))
})
print(do.call(rbind, res), right = FALSE)
## variance SE mapping: hand delta-check on VVI, EEI by observed information refit on variance param
f <- multilpa(d, vars[1:3], id = "student", n_profiles = 2, n_group_classes = 1, model = "EEI", n_starts = 2, seed = 1)
print(get_results(f, "profiles")); print(subset(as.data.frame(parameter_inference(f)), parameter %in% c("variance","covariance")))
## Table vs plot on a two-level fit
f2 <- multilpa(course_engagement, vars, id = "student", n_profiles = 3, n_group_classes = 2, n_starts = 3, seed = 1)
tp <- get_results(f2, "profiles"); ts <- get_results(f2, "profiles", scale = "standardized")
hm <- plot_heatmap(f2); hd <- ld(hm, 1)
cat("heatmap fill scale & data vs table: \n")
hm_data <- hm$data; ind <- as.character(hm_data$indicator)
chk <- merge(hm_data[c("profile","position","mean")], data.frame(profile = ts$profile, position = match(ts$indicator, vars), tmean = ts$mean))
print(all.equal(chk$mean, chk$tmean))
# standardized by hand
gd <- get_results(f2, "data")
hand <- (tp$mean - colMeans(gd[vars])[tp$indicator]) / apply(gd[vars], 2, sd)[tp$indicator]
print(all.equal(unname(hand), ts$mean))
print(all.equal(ts$mean_standard_error, unname(tp$mean_standard_error / apply(gd[vars], 2, sd)[tp$indicator])))
## bar / profile plots vs table, labels vs counts
pb <- plot_bars(f2); bd <- pb$data
print(all.equal(bd$mean[order(bd$profile, bd$position)], tp$mean))
print(all.equal(((bd$upper-bd$lower)/(2*qnorm(.975)))[order(bd$profile, bd$position)], tp$mean_standard_error))
print(levels(bd$profile_label)); print(subset(get_results(f2,"counts"), level == "individuals"))
# sizes, entropy, posterior, avepp list layer vs tables
ps <- plot_sizes(f2); print(class(ps)); print(names(ps))
print(ps$individuals$data[c("profile","count","share","bar_label")]); print(ps$groups$data[c("profile","count","share","bar_label")])
pe <- plot_entropy(f2); print(pe$individuals$labels$subtitle); print(get_entropy(f2))
cat("mean case entropy vs 1-rel:", mean(pe$individuals$data$value), 1 - get_entropy(f2)$relative_entropy[1], "\n")
pp <- plot_posteriors(f2); ppd <- pp$individuals$data; print(tapply(ppd$value, ppd$profile, mean)); print(subset(get_classification(f2), level == "individuals")[c("class","n_modal","average_posterior")]); 
pa <- plot_avepp(f2); print(pa$individuals$data[c("row","column","average_posterior")]); print(head(get_avepp(f2), 9))
print(pa$individuals$data$assigned_label[1:3])
## prevalence plot vs table
pv <- plot_prevalence(f2); print(pv$data[c("position","class","mean","profile_label")]); print(get_results(f2,"profile_probabilities")); print(pv$scales$get_scales("x")$labels)
## transitions plot vs table
mv <- lta(course_engagement, vars, "student", n_profiles = 3, n_group_classes = 2, time = "sequence", n_starts = 3, seed = 1)
pt <- plot_transitions(mv); ptd <- pt$data
tt <- get_results(mv, "transitions"); m <- merge(ptd[c("from","to","class","probability")], tt, by.x = c("from","to","class"), by.y = c("from","to","group_class"))
print(all.equal(m$probability.x, m$probability.y))
## Non-profile families: wild plots vs table
reg <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2, n_starts = 3, seed = 1)
print(get_results(reg, "classes"))
s <- plot_sizes(reg); print(class(s)); print(s$data[c("class","count","share")])
print(get_classification(reg)); print(get_entropy(reg))
e <- plot_entropy(reg); print(e$labels$subtitle); print(e$labels$x); cat("1-mean:", 1-mean(e$data$value), "\n")
po <- plot_posteriors(reg); print(tapply(po$data$value, po$data$class, mean))
av <- plot_avepp(reg); print(av$data[c("row","column","average_posterior")]); print(get_avepp(reg))
gr <- mixture_regression(score ~ wave, subset(growth_scores, student <= 80), n_classes = 2, id = "student", class_level = "group", random = "wave", n_starts = 2, seed = 1)
print(get_classification(gr)); print(plot_sizes(gr)$data[c("class","count","share")])
print(get_results(gr,"classes"))
print(get_entropy(gr)); ge <- plot_entropy(gr); print(ge$labels$subtitle)
