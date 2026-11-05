# Independent synthetic probes for the 2026-11-05 review.
pkgload::load_all(".", quiet = TRUE)
options(width = 120)
warned <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, paste(class(w)[1L], conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  attr(value, "review_warnings") <- warnings
  value
}
section <- function(title) cat("\n", title, "\n", sep = "")

section("E1 partial singleton centring")
set.seed(110501)
sizes <- c(rep(1L, 10), rep(6L, 65))
id <- rep(seq_along(sizes), sizes)
state <- sample(1:2, length(id), TRUE)
offset <- rep(rnorm(length(sizes), sd = 1.4), sizes)
dd <- data.frame(id = id, a = offset + c(-1, 1)[state] + rnorm(length(id), sd = .65),
                 b = offset + c(1, -1)[state] + rnorm(length(id), sd = .65),
                 c = rnorm(length(id), c(-.4, .4)[state]))
str(dd)
fits <- lapply(2:3, function(k) warned(multilpa(dd, c("a", "b", "c"), "id", k, 1,
  centering = "person", n_starts = 4, seed = 193, tol = 1e-10, cores = 1)))
print(data.frame(k = 2:3, bic = vapply(fits, function(f) f$bic, 0),
                 boundary = vapply(fits, function(f) f$boundary, FALSE)))
print(fits[[2]]$means)
print(fits[[2]]$variances)
cat("All-zero rows:", sum(rowSums(abs(fits[[2]]$indicator_data)) == 0), "\n")
print(attr(fits[[2]], "review_warnings"))
trimmed <- dd[dd$id > 10L, ]
trimmed_fits <- lapply(2:3, function(k) warned(multilpa(trimmed, c("a", "b", "c"), "id", k, 1,
  centering = "person", n_starts = 4, seed = 193, tol = 1e-10, cores = 1)))
print(data.frame(k = 2:3, trimmed_bic = vapply(trimmed_fits, function(f) f$bic, 0)))
stopifnot(sum(rowSums(abs(fits[[2]]$indicator_data)) == 0) == 10L)

section("E2 weighted grand transformations against duplicated rows")
set.seed(110502)
n <- 320L
state <- sample(0:1, n, TRUE)
ww <- data.frame(a = rnorm(n, 3 * state), b = rnorm(n, -2 * state, 1 + state),
                 c = rnorm(n, state), w = ifelse(state == 0, 3L, 1L))
duplicate <- ww[rep(seq_len(n), ww$w), ]
compare_weighted <- function(transform) {
  fw <- warned(lpa(ww, c("a", "b", "c"), 2, model = "EII", weights = "w",
    centering = transform, scaling = transform, n_starts = 3, seed = 14, tol = 1e-10, cores = 1))
  fd <- warned(lpa(duplicate, c("a", "b", "c"), 2, model = "EII",
    centering = transform, scaling = transform, n_starts = 3, seed = 14, tol = 1e-10, cores = 1))
  raw <- function(f) {
    value <- f$means[order(f$means[, "a"]), , drop = FALSE]
    if (!is.null(f$scaling_factors)) value <- sweep(value, 2, f$scaling_factors[1, ], "*")
    if (!is.null(f$centering_offsets)) value <- sweep(value, 2, f$centering_offsets[1, ], "+")
    value
  }
  if (transform == "grand") {
    print(fw$centering_offsets)
    print(fd$centering_offsets)
    print(fw$scaling_factors)
    print(fd$scaling_factors)
  }
  data.frame(transform = transform, raw_mean_difference = max(abs(raw(fw) - raw(fd))),
    mean_difference = max(abs(fw$means[order(fw$means[, 1]), ] - fd$means[order(fd$means[, 1]), ])))
}
print(do.call(rbind, lapply(c("none", "grand"), compare_weighted)))

section("E3 membership unused level")
set.seed(110503)
person <- seq_len(70L)
rr <- data.frame(id = rep(person, each = 4L), t = rep(0:3, length(person)))
state <- rep(sample(0:1, length(person), TRUE), each = 4)
rr$grp <- factor(rep(ifelse(person %% 2 == 0, "a", "b"), each = 4), levels = c("a", "b", "unused"))
rr$y <- 3 * state + (.8 - 1.3 * state) * rr$t + rnorm(nrow(rr), sd = .65)
fm <- warned(mixture_regression(y ~ t, rr, 2, id = "id", class_level = "group",
  membership = ~ grp, n_starts = 3, seed = 53, tol = 1e-10))
rr_drop <- droplevels(rr)
fd <- warned(mixture_regression(y ~ t, rr_drop, 2, id = "id", class_level = "group",
  membership = ~ grp, n_starts = 3, seed = 53, tol = 1e-10))
print(data.frame(case = c("unused retained", "unused dropped"),
                 k = c(fm$n_parameters, fd$n_parameters),
                 LL = c(fm$log_likelihood, fd$log_likelihood),
                 bic = c(get_results(fm, "fit")$bic, get_results(fd, "fit")$bic)))
print(colnames(fm$spec$w)); print(attr(fm, "review_warnings"))
stopifnot(fm$n_parameters == fd$n_parameters + 1,
          abs(fm$log_likelihood - fd$log_likelihood) < 1e-5)

section("E4 reserved id class overwrites ids")
reserved <- rr_drop
reserved$class <- reserved$id + 1000L
f_reserved <- warned(mixture_regression(y ~ t, reserved, 2, id = "class", class_level = "group",
  n_starts = 2, seed = 41))
group_table <- as.data.frame(get_results(f_reserved, "groups"))
print(head(group_table)); print(names(group_table)); print(attr(f_reserved, "review_warnings"))
cat("Input IDs range:", range(reserved$class), "reported class range:", range(group_table$class), "\n")
stopifnot(all(group_table$class %in% 1:2), !any(group_table$class %in% reserved$class))

section("E5 growth two-occasion exact covariance ridge")
set.seed(110504)
growth <- data.frame(id = rep(1:120, each = 2), t = rep(0:1, 120))
growth$y <- 1 + .7 * growth$t + rep(rnorm(120), each = 2) + rnorm(nrow(growth), sd = .7)
fg <- warned(mixture_regression(y ~ t, growth, 1, id = "id", class_level = "group",
  random = ~ 1 + t, n_starts = 1, seed = 16, tol = 1e-8))
params <- fg$params
ridge <- do.call(rbind, lapply(c(-.05, 0, .025), function(delta) {
  p <- params
  G <- p$random_covariance[[1]]
  G[1, 1] <- G[1, 1] - delta
  G[1, 2] <- G[2, 1] <- G[1, 2] + delta
  G[2, 2] <- G[2, 2] - 2 * delta
  p$random_covariance <- list(G)
  p$sigma2 <- p$sigma2 + delta
  data.frame(delta = delta, min_eigen = min(eigen(G, symmetric = TRUE)$values),
    residual_variance = p$sigma2, LL = .growth_expectation(fg$spec, fg$stats, p)$log_likelihood)
}))
print(ridge)
cat("Converged:", fg$converged, "n_parameters:", fg$n_parameters, "random boundary:", fg$random_boundary, "\n")
print(attr(fg, "review_warnings"))
stopifnot(all(ridge$min_eigen > 0), diff(range(ridge$LL)) < 1e-8)

section("E6 occasion transitions bootstrap loses unique longest person")
set.seed(110505)
lt <- data.frame(id = rep(1:50, each = 3), time = rep(1:3, 50))
state <- sample(0:1, nrow(lt), TRUE)
lt$a <- rnorm(nrow(lt), 4 * state, .7)
lt$b <- rnorm(nrow(lt), -3 * state, .7)
lt <- rbind(lt, data.frame(id = 1L, time = 4:5, a = c(.2, 3.7), b = c(.1, -2.8)))
fl <- warned(lta(lt, c("a", "b"), "id", 2, "time", transitions = "occasion",
  n_starts = 2, seed = 52, tol = 1e-8, cores = 1))
boot <- tryCatch(warned(parameter_inference(fl, data = lt, method = "bootstrap", iter = 15,
  n_starts = 1, max_iter = 1000, seed = 48)), error = function(e) e)
cat("Bootstrap result class:", class(boot), "\n")
if (inherits(boot, "error")) cat("Bootstrap error:", conditionMessage(boot), "\n") else print(attributes(boot)[c("n_valid", "messages")])
saveRDS(list(singletons = fits, trimmed = trimmed_fits, membership = fm, membership_dropped = fd,
  reserved = f_reserved, growth = fg, growth_ridge = ridge, transitions = fl, bootstrap_result = boot),
  "validation/review-2026-11-05/engines/results.rds")
