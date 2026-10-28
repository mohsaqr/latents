pkgload::load_all(quiet = TRUE)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
capture <- function(expr) tryCatch(expr, error = function(e) list(error = conditionMessage(e), class = class(e)))
results <- list()
set.seed(101)
d <- data.frame(g = rep(1:20, each = 8), t = rep(1:8, 20), y = rnorm(160), z = rnorm(160))
grid <- quiet(enumerate_classes(d, c("y", "z"), "g", n_profiles = 2:3,
                               n_group_classes = 1, model = "VVI", n_starts = 2, seed = 1))
cmp <- compare_models(grid)
small <- head(cmp, 1)
results$comparison_subset <- list(present = small$model,
                                 classification = unique(get_results(small, "classification")$model),
                                 attributes_models = names(attr(small, "models")))
results$comparison_all_mixture_regression <- capture({
 f1 <- quiet(mixture_regression(y ~ z, d, n_classes = 1, seed = 1))
 f2 <- quiet(mixture_regression(y ~ z, d, n_classes = 2, seed = 1))
 get_results(compare_models(a = f1, b = f2), "all")
})
crossed <- quiet(enumerate_classes(d, c("y", "z"), "g", n_profiles = 2:3,
                                  n_group_classes = 1, model = c("EEI", "VVI"),
                                  n_starts = 2, seed = 1))
p <- plot(crossed, what = "smallest")
results$smallest_crossed <- list(candidates = as.data.frame(crossed)[c("candidate", "model")],
                                plot = p$data,
                                built_groups = ggplot2::ggplot_build(p)$data[[2]][c("x", "y", "group", "PANEL")])
trans <- quiet(enumerate_classes(d, c("y", "z"), "g", time = "t", n_profiles = 2:3,
                                n_group_classes = 1, n_starts = 1, seed = 1))
results$direct_transition <- capture(compare_models(a = get_candidate(trans, 2), b = get_candidate(trans, 3)))
results$transition_all <- capture(get_results(trans, "all"))
results$transition_null <- capture(get_results(trans, NULL))
results$transition_grid <- capture(compare_models(trans))
general <- quiet(enumerate_classes(d, c("y", "z"), "g", time = "t", n_profiles = 2:3,
                                   n_group_classes = 1, transitions = "occasion", n_starts = 1, seed = 1))
results$general_grid <- capture(compare_models(general))
fam <- quiet(enumerate_classes(d, c("y", "z"), "g", family = "additive",
                              n_group_classes = 1:2, n_starts = 1, seed = 1))
results$family_all <- capture(get_results(fam, "all"))
results$family_null <- capture(get_results(fam, NULL))
# Identical centered values with different raw observations.
d2 <- transform(d, y = y + 100, z = z + 100)
f1 <- quiet(multilpa(d, c("y", "z"), "g", n_profiles = 2, n_group_classes = 1,
                    centering = "grand", n_starts = 1, seed = 1))
f2 <- quiet(multilpa(d2, c("y", "z"), "g", n_profiles = 3, n_group_classes = 1,
                    centering = "grand", n_starts = 1, seed = 1))
results$centered_different_data <- capture(compare_models(a = f1, b = f2))
# Different categorical values mapped to identical integer codes.
dcat <- transform(d, a = ifelse(y > 0, "low", "high"), b = ifelse(z > 0, "low", "high"))
dcat2 <- transform(dcat, a = ifelse(a == "low", "yes", "no"), b = ifelse(b == "low", "yes", "no"))
c1 <- quiet(lca(dcat, c("a", "b"), n_classes = 2, n_starts = 1, seed = 1))
c2 <- quiet(lca(dcat2, c("a", "b"), n_classes = 2, n_starts = 1, seed = 1))
results$categorical_different_data <- capture(compare_models(a = c1, b = c2))
saveRDS(results, "validation/audit-2026-10-13/recent-api/probes.rds")
str(results, max.level = 3)
