# Which kinds of start reach which maximum on heavy_02.
suppressMessages(pkgload::load_all(".", quiet = TRUE))
invisible(lapply(c("helper-montecarlo.R", "helper-mplus-classification.R"), \(file) {
  sys.source(file.path("equivalence", file), envir = globalenv())
}))
root <- file.path("equivalence", "features", "direct-effects")
sys.source(file.path(root, "helpers.R"), envir = globalenv())
case <- de_cases(root)$heavy_02
vars <- case$spec$vars
fit_with <- \(search, n) suppressWarnings(lpa(case$data, vars, 2, model = "VVI", n_starts = n,
  seed = 1, tol = 1e-10, max_iter = 5000, direct_effects = "y3 ~ z", search = search))
basin <- \(ll) cut(ll, c(-Inf, -3106.99, -3106.98, -3098.49, -3098.48, Inf),
                   c("below low", "low -3106.986", "between", "high -3098.487", "above"))
km <- fit_with("kmeans", 200)
st <- get_results(km, "starts")
str(st)
st$basin <- basin(st$log_likelihood)
st$kind <- vapply(st$start, \(i) .latents_start_kind_label(i, 200L), character(1))
print(with(st, table(kind, basin, useNA = "ifany")))
print(head(st, 3))
em <- fit_with("emEM", 200)
se <- get_results(em, "starts")
str(se)
se$basin <- basin(se$log_likelihood)
se$kind <- vapply(se$start, \(i) .latents_start_kind_label(i, 200L), character(1))
print(with(subset(se, stage == "final"), table(kind, basin, useNA = "ifany")))
