.libPaths(c("/private/tmp/latents-review-20261022/rlib", .libPaths()))
library(latents)
set.seed(14)
data <- data.frame(y1 = rnorm(70, rep(c(-2, 2), each = 35)),
                   y2 = rnorm(70, rep(c(-1, 1), each = 35)),
                   y3 = rnorm(70), group = rep(1:14, each = 5))
cases <- expand.grid(model = latents:::.multilpa_structures(),
                      acceleration = c("none", "squarem"), groups = c(1L, 2L),
                      stringsAsFactors = FALSE)
results <- lapply(seq_len(nrow(cases)), function(i) {
  args <- c(list(data = data, vars = c("y1", "y2", "y3"), id = "group",
                n_profiles = 2, n_starts = 2, max_iter = 40, seed = 3,
                model = cases$model[i], acceleration = cases$acceleration[i],
                n_group_classes = cases$groups[i]))
  fit <- function(search) suppressWarnings(do.call(multilpa, c(args,
    list(search = search), if (search == "emEM") list(n_screen = 2L))))
  a <- tryCatch(fit("kmeans"), error = identity)
  b <- tryCatch(fit("emEM"), error = identity)
  list(case = cases[i, ], error = c(inherits(a, "error"), inherits(b, "error")),
       same = identical(a$log_likelihood, b$log_likelihood) &&
         identical(a$means, b$means) && identical(a$covariances, b$covariances) &&
         identical(a$history, b$history))
})
saveRDS(results, "/private/tmp/latents-review-20261022/probe-emem-identity.rds")
cat("Cases", length(results), "nonidentical", sum(!vapply(results, `[[`, logical(1), "same")), "\n")
print(Filter(function(x) !x$same || any(x$error), results))
