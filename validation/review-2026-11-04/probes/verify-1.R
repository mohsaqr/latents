# Parent's independent reproduction of the MAJOR findings, on new seeds.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL)
q <- function(e) { w <- character(); v <- withCallingHandlers(tryCatch(e, error = function(x) structure(list(msg = conditionMessage(x), cls = class(x)[1]), class = "failed")),
  warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") }); attr(v, "w") <- unique(w); v }
sec <- function(s) cat("\n=====", s, "\n")
one <- function(label, code) { sec(label); r <- tryCatch(code, error = function(e) cat("SCRIPT ERROR:", conditionMessage(e), "\n")); invisible(r) }

one("C1 person centring with singleton groups", {
  for (s in c(101, 102, 103)) { set.seed(s)
    J <- 80; sizes <- c(rep(1L, 12), sample(5:9, J - 12, TRUE)); g <- rep(seq_len(J), sizes); n <- length(g)
    z <- rbinom(n, 1, .5); off <- rep(rnorm(J, 0, 1.5), sizes)
    d <- data.frame(g = g, a = off + rnorm(n, c(-1, 1)[z + 1], .7), b = off + rnorm(n, c(1, -1)[z + 1], .7), c = rnorm(n, c(0, .5)[z + 1], 1))
    d2 <- subset(d, g %in% which(sizes > 1)); v <- c("a", "b", "c")
    fit <- function(D, k) q(multilpa(D, v, "g", n_profiles = k, n_group_classes = 1, centering = "person", n_starts = 5, seed = 1))
    a2 <- fit(d, 2); a3 <- fit(d, 3); b2 <- fit(d2, 2); b3 <- fit(d2, 3)
    cat(sprintf("seed %d: with singletons BIC k2 %.1f k3 %.1f | without BIC k2 %.1f k3 %.1f | k3 min variance %.2g, warnings: %s\n", s, a2$bic, a3$bic, b2$bic, b3$bic, min(a3$variances), paste(attr(a3, "w"), collapse = ","))) }
})

one("C2 weights with grand centring/scaling vs duplicated rows", {
  set.seed(201); n <- 400; z <- rbinom(n, 1, .5)
  d <- data.frame(a = rnorm(n, z * 3, 1), b = rnorm(n, -z * 2, c(1, 2)[z + 1]), c = rnorm(n, z, 1)); d$w <- ifelse(z == 1, 1L, 3L)
  dup <- d[rep(seq_len(n), d$w), ]; dup$w <- NULL; rownames(dup) <- NULL; v <- c("a", "b", "c")
  for (cfg in list(c("VVI", "none", "none"), c("EII", "none", "none"), c("VVI", "grand", "none"), c("VVI", "grand", "grand"), c("EII", "grand", "grand"))) {
    fw <- q(lpa(d, v, 2, model = cfg[1], weights = "w", centering = cfg[2], scaling = cfg[3], n_starts = 5, seed = 1, tol = 1e-10))
    fd <- q(lpa(dup, v, 2, model = cfg[1], centering = cfg[2], scaling = cfg[3], n_starts = 5, seed = 1, tol = 1e-10))
    o1 <- order(fw$means[, 1]); o2 <- order(fd$means[, 1])
    cat(sprintf("%s c=%s s=%s: max|mean diff| %.5f  shares diff %.5f\n", cfg[1], cfg[2], cfg[3], max(abs(fw$means[o1, ] - fd$means[o2, ])), max(abs(sort(c(fw$profile_probabilities)) - sort(c(fd$profile_probabilities)))))) }
})

one("I1 two-step LCA: held response logits in the information", {
  set.seed(301); n <- 500; z <- rnorm(n); k <- rbinom(n, 1, plogis(.2 + .8 * z)) + 1L; p <- c(.75, .3)
  d <- data.frame(g = seq_len(n), z = z); v <- paste0("u", 1:5); for (x in v) d[[x]] <- factor(ifelse(runif(n) < p[k], "yes", "no"))
  base <- multilpa(d, v, "g", n_profiles = 2, n_group_classes = 1, categorical = v, n_starts = 5, seed = 1, tol = 1e-10, max_iter = 5000)
  two <- multilpa(d, v, "g", n_profiles = 2, n_group_classes = 1, categorical = v, profile_covariates = "z", start = starting_values(base, what = "measurement"), fixed = "measurement", n_starts = 1, seed = 1, tol = 1e-10, max_iter = 5000)
  rp <- two$response_probabilities
  ld <- vapply(1:2, function(j) rowSums(vapply(v, function(x) log(rp[[x]][j, as.character(d[[x]])]), numeric(n))), numeric(n))
  ll <- function(b) { P <- plogis(b[1] + b[2] * z); sum(log(P * exp(ld[, 1]) + (1 - P) * exp(ld[, 2]))) }
  h <- function(f, x, e = 1e-4) { m <- length(x); H <- matrix(0, m, m); for (i in seq_len(m)) for (j in seq_len(m)) { a <- b <- c0 <- d0 <- x; a[i] <- a[i] + e; a[j] <- a[j] + e; b[i] <- b[i] + e; b[j] <- b[j] - e; c0[i] <- c0[i] - e; c0[j] <- c0[j] + e; d0[i] <- d0[i] - e; d0[j] <- d0[j] - e; H[i, j] <- (f(a) - f(b) - f(c0) + f(d0)) / (4 * e^2) }; H }
  cf <- as.vector(two$profile_coefficients)
  cat("n_parameters", two$n_parameters, " package SE", signif(parameter_inference(two)$standard_error, 5), " conditional SE by own finite differences", signif(sqrt(diag(solve(-h(ll, cf)))), 5), "\n")
  cat("vcov unconstrained dim", dim(vcov(two, scale = "unconstrained")), "\n")
})

one("I2 OPG SE of a mixing probability under grand transforms", {
  set.seed(401); G <- 50; m <- 6; h <- rep(rbinom(G, 1, .5), each = m); k <- rbinom(G * m, 1, plogis(-.5 + 1.5 * h)) + 1L
  d <- data.frame(g = rep(seq_len(G), each = m), y1 = rnorm(G * m, c(-1.5, 1.5)[k]), y2 = rnorm(G * m, c(1, -1)[k]), y3 = rnorm(G * m, c(0, 1)[k], .7))
  a <- list(data = d, vars = c("y1", "y2", "y3"), id = "g", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1, tol = 1e-10, max_iter = 5000)
  f0 <- do.call(multilpa, a); fc <- do.call(multilpa, c(a, list(centering = "grand"))); fs <- do.call(multilpa, c(a, list(centering = "grand", scaling = "grand")))
  for (ty in c("observed", "robust", "opg")) { se <- function(f) { p <- parameter_inference(f, vcov_type = ty); signif(p$standard_error[p$parameter == "probability"], 5) }
    cat(ty, "none:", se(f0), "| centred:", se(fc), "| both:", se(fs), "\n") }
})

one("S1 boundary flag on LCA at the floor", {
  set.seed(501); n <- 300; k <- rbinom(n, 1, .5) + 1
  d <- data.frame(a = ifelse(k == 1, "x", ifelse(runif(n) < .5, "x", "y")), b = ifelse(k == 1, "u", ifelse(runif(n) < .3, "u", "v")), c = ifelse(runif(n) < c(.8, .2)[k], "r", "s"), d = ifelse(runif(n) < c(.9, .3)[k], "p", "q"), e = ifelse(runif(n) < c(.7, .2)[k], "m", "n"))
  l <- q(lca(d, letters[1:5], n_classes = 2, n_starts = 6, seed = 1))
  cat("min prob", min(unlist(l$response_probabilities)), " floor", .multilpa_probability_floor(l), " fit$boundary", l$boundary, "\n")
  ch <- get_checks(l); print(subset(ch, grepl("bound", check, ignore.case = TRUE))[, intersect(names(ch), c("check", "value", "status"))], row.names = FALSE)
  w <- q(parameter_inference(l)); cat("parameter_inference:", if (inherits(w, "failed")) paste("REFUSED", w$cls) else "ran", "\n")
})

one("S2 bootstrap_lrt with a categorical joint pair", {
  set.seed(601); ng <- 40; np <- 5; g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each = np); k <- 1 + rbinom(length(id), 1, ifelse(g[id] == 1, .8, .2))
  d <- data.frame(id = id, y1 = rnorm(length(id), c(-1.5, 1.5)[k])); d$a <- ifelse(runif(nrow(d)) < c(.2, .8)[k], "x", "y"); d$b <- ifelse(runif(nrow(d)) < c(.7, .3)[k], "u", "v")
  mk <- function(H) q(multilpa(d, c("a", "b", "y1"), "id", n_profiles = 2, n_group_classes = H, categorical = c("a", "b"), local_dependence = "a:b", model = "VVI", n_starts = 2, seed = 1))
  r <- q(bootstrap_lrt(mk(1), mk(2), iter = 3, seed = 1))
  if (inherits(r, "failed")) cat("FAILED", r$cls, r$msg, "\n") else { print(get_results(r, "test")[, intersect(names(get_results(r, "test")), c("statistic", "lr", "p_value", "n_valid", "iter"))]); cat("warnings:", attr(r, "w"), "\n"); print(unique(get_results(r, "replicates")$error)) }
})

one("R1 random intercept and slope on two occasions", {
  set.seed(701); pp <- do.call(rbind, lapply(1:150, function(i) { b <- rnorm(1); data.frame(id = i, t = 0:1, y = b + .5 * (0:1) + rnorm(2, sd = .7)) }))
  f <- q(mixture_regression(y ~ t, pp, 1, id = "id", class_level = "group", random = ~ 1 + t, seed = 1))
  p <- f$params; G <- p$random_covariance[[1]]; s2 <- p$sigma2
  at <- function(dl) { a <- p; g <- G; g[1, 1] <- g[1, 1] - dl; g[1, 2] <- g[2, 1] <- g[1, 2] + dl; g[2, 2] <- g[2, 2] - 2 * dl; a$random_covariance <- list(g); a$sigma2 <- s2 + dl; .growth_expectation(f$spec, f$stats, a)$log_likelihood }
  # own marginal likelihood, without package code, at d = 0 and d = 0.03
  own <- function(dl) { g <- G; g[1, 1] <- g[1, 1] - dl; g[1, 2] <- g[2, 1] <- g[1, 2] + dl; g[2, 2] <- g[2, 2] - 2 * dl; Z <- cbind(1, 0:1); V <- Z %*% g %*% t(Z) + diag(s2 + dl, 2); mu <- drop(cbind(1, 0:1) %*% as.vector(p$beta %||% p$coefficients)); sum(vapply(split(pp$y, pp$id), function(y) { e <- y - mu; -log(2 * pi) - .5 * determinant(V)$modulus[1] - .5 * sum(e * solve(V, e)) }, 0)) }
  cat("LL at d=-0.1, 0, 0.03 (package):", format(c(at(-.1), at(0), at(.03)), digits = 10), "\n")
  r <- tryCatch(format(c(own(0), own(.03)), digits = 10), error = function(e) paste("own failed", conditionMessage(e))); cat("own LL at d=0, .03:", r, "\n")
  cat("converged", f$converged, "random_boundary", f$random_boundary, "n_parameters", f$n_parameters, "warnings:", attr(f, "w"), "\n")
})

one("R2 unused factor level in membership", {
  set.seed(801); d <- do.call(rbind, lapply(1:60, function(i) { gr <- sample(c("a", "b"), 1); k <- rbinom(1, 1, ifelse(gr == "a", .7, .3)); data.frame(id = i, grp = gr, t = 0:3, y = c(0, 4)[k + 1] + c(.5, -.5)[k + 1] * (0:3) + rnorm(4, sd = .8)) }))
  d1 <- d; d1$grp <- factor(d1$grp, levels = c("a", "b", "zz")); d2 <- d; d2$grp <- factor(d2$grp)
  m <- function(D) q(mixture_regression(y ~ t, D, 2, id = "id", class_level = "group", membership = ~ grp, n_starts = 3, seed = 1))
  a <- m(d1); b <- m(d2); cat(sprintf("unused level: k %d LL %.4f BIC %.3f warnings %s\ndropped:      k %d LL %.4f BIC %.3f warnings %s\n", a$n_parameters, a$log_likelihood, a$bic, paste(attr(a, "w"), collapse = ","), b$n_parameters, b$log_likelihood, b$bic, paste(attr(b, "w"), collapse = ",")))
  cat("SE all NA with unused level:", all(is.na(get_results(a, "membership")$standard_error)), "; without:", all(is.na(get_results(b, "membership")$standard_error)), "\n")
})

one("R3 column named class as id", {
  set.seed(901); d <- do.call(rbind, lapply(1:40, function(i) { k <- rbinom(1, 1, .5); data.frame(class = 500 + i, t = 0:3, y = c(0, 4)[k + 1] + rnorm(4, sd = .8)) }))
  f <- q(mixture_regression(y ~ t, d, 2, id = "class", class_level = "group", n_starts = 2, seed = 1)); g <- get_results(f, "groups"); print(names(g)); print(head(g, 3), row.names = FALSE)
})

one("M1 compare_models joint against staged", {
  v <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
  f0 <- q(multilpa(course_engagement, v, id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1)); st <- q(fit_staged(course_engagement, v, id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 3, seed = 1))
  print(as.data.frame(compare_models(joint = f0, staged = st))[c("model", "n_parameters", "log_likelihood", "bic_groups", "delta_bic", "bic_weight")], row.names = FALSE)
  cat("staged n_parameters_with_measurement", st$n_parameters_with_measurement, " BIC with it", -2 * st$log_likelihood + st$n_parameters_with_measurement * log(st$n_groups), "\n")
})

one("M2 p_adjusted in the two residual tables", {
  f <- q(lca(student_esm, c("time_with_friends", "on_social_media", "sports", "listened_music", "walking"), n_classes = 3, n_starts = 3, seed = 1))
  a <- get_results(f, "residuals"); b <- get_residuals(f); cat("get_results p_adjusted == p_value:", isTRUE(all.equal(a$p_adjusted, a$p_value)), "; get_residuals:", isTRUE(all.equal(b$p_adjusted, b$p_value)), "\n")
  cat("pairs below .01: get_results", sum(a$p_adjusted < .01, na.rm = TRUE), " get_residuals", sum(b$p_adjusted < .01, na.rm = TRUE), "\n")
})
