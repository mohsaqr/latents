suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL)
q <- function(e) { w <- character(); v <- withCallingHandlers(tryCatch(e, error = function(x) structure(list(msg = conditionMessage(x), cls = class(x)[1]), class = "failed")),
  warning = function(c) { w <<- c(w, class(c)[1]); invokeRestart("muffleWarning") }); attr(v, "w") <- unique(w); v }
cat("===== R2c\n")
set.seed(801); d <- do.call(rbind, lapply(1:60, function(i) { gr <- sample(c("a", "b"), 1); k <- rbinom(1, 1, ifelse(gr == "a", .7, .3)); data.frame(id = i, grp = gr, t = 0:3, y = c(0, 4)[k + 1] + c(.5, -.5)[k + 1] * (0:3) + rnorm(4, sd = .8)) }))
d1 <- d; d1$grp <- factor(d1$grp, levels = c("a", "b", "zz")); d2 <- d; d2$grp <- factor(d2$grp)
m <- function(D) q(mixture_regression(y ~ t, D, 2, id = "id", class_level = "group", membership = ~ grp, n_starts = 3, seed = 1))
a <- m(d1); b <- m(d2); fa <- as.data.frame(get_results(a, "fit")); fb <- as.data.frame(get_results(b, "fit"))
print(rbind(unused = fa[, c("n_parameters", "log_likelihood", "bic")], dropped = fb[, c("n_parameters", "log_likelihood", "bic")])); cat("warnings unused:", attr(a, "w"), "| dropped:", attr(b, "w"), "\n")
cat("===== T3b\n")
set.seed(1501); P <- matrix(c(.8, .2, .3, .7), 2, byrow = TRUE); path <- function(i, TT) { s <- integer(TT); s[1] <- sample(2, 1, prob = c(.6, .4)); for (t in 2:TT) s[t] <- sample(2, 1, prob = P[s[t - 1], ]); data.frame(id = i, time = seq_len(TT), y1 = rnorm(TT, c(0, 3)[s]), y2 = rnorm(TT, c(1, -2)[s], .8)) }
dd <- do.call(rbind, c(lapply(1:60, path, TT = 3), list(path(61, 5))))
o <- q(lta(dd, c("y1", "y2"), "id", 2, "time", transitions = "occasion", n_starts = 2, seed = 1)); h <- q(lta(dd, c("y1", "y2"), "id", 2, "time", n_starts = 2, seed = 1))
b1 <- q(parameter_inference(o, data = dd, method = "bootstrap", iter = 30, n_starts = 1, seed = 1)); if (inherits(b1, "failed")) cat("occasion FAILED:", b1$cls, "|", b1$msg, "\n") else cat("occasion completed", nrow(b1), "\n")
b2 <- q(parameter_inference(h, data = dd, method = "bootstrap", iter = 30, n_starts = 1, seed = 1)); if (inherits(b2, "failed")) cat("homogeneous FAILED:", b2$cls, "|", b2$msg, "\n") else cat("homogeneous completed", nrow(b2), "\n")
cat("===== drop_empty, every criterion, one empty row inside a group and one group made only of empty rows\n")
set.seed(1601); G <- 30; m <- 4; g <- rep(seq_len(G), each = m); k <- rbinom(G * m, 1, .5); v <- c("a", "b", "c")
dm <- data.frame(g = g, a = rnorm(G * m, c(-2, 2)[k + 1]), b = rnorm(G * m, c(1, -1)[k + 1]), c = rnorm(G * m)); dm[1, v] <- NA; dm[dm$g == 2, v] <- NA
t1 <- multilpa(dm, v, "g", n_profiles = 2, n_group_classes = 2, missing = "fiml", drop_empty = TRUE, n_starts = 3, seed = 1); t2 <- multilpa(dm, v, "g", n_profiles = 2, n_group_classes = 2, missing = "fiml", drop_empty = FALSE, n_starts = 3, seed = 1)
i1 <- get_results(t1, "information_criteria"); i2 <- get_results(t2, "information_criteria"); print(data.frame(criterion = names(i1), true_ = round(unlist(i1), 4), false_ = round(unlist(i2), 4), differs = abs(unlist(i1) - unlist(i2)) > 1e-8), row.names = FALSE)
cat("logLik equal:", isTRUE(all.equal(t1$log_likelihood, t2$log_likelihood)), "\n"); print(get_results(t1, "entropy")); print(get_results(t2, "entropy"))
