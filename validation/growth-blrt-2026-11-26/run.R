# Runs the eval = FALSE `blrt` chunk of vignettes/articles/growth-mixture.Rmd
# once, so the p-values the article reports come from a recorded run.
pkgload::load_all(".", quiet = TRUE)
options(width = 120)
timing <- system.time(
  grid <- enumerate_regressions(score ~ wave, growth_scores, n_classes = 1:4,
                                id = "student", class_level = "group",
                                random = "wave", random_covariance = "equal",
                                bootstrap = 99, seed = 1)
)
print(timing)
print(grid)
print(as.data.frame(grid))
print(sessionInfo())
