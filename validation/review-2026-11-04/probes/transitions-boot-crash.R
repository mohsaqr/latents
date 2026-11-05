suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
sim <- function(n, T, seed, K2 = TRUE) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 3)[s]), y2 = rnorm(T, c(1, -2)[s], .8))
  }))
}
vars <- c("y1","y2")
d3 <- sim(60, 3, 6); extra <- transform(subset(sim(1, 5, 7), time > 3), id = 1L); d3 <- rbind(d3, extra)
ids <- table(d3$id); cat("persons", length(ids), " lengths", names(table(ids)), table(ids), "\n")
o <- suppressWarnings(lta(d3, vars, "id", 2, "time", transitions = "occasion", n_starts = 2, seed = 1))
cat("n_occasions", o$n_occasions, " terms", dimnames(o$transition_coefficients)[[1]], "\n")
# What each replicate looks like
set.seed(1)
rows <- split(seq_len(nrow(d3)), o$group_index)
res <- vapply(1:12, function(i) {
  draw <- sample.int(o$n_groups, o$n_groups, TRUE)
  rs <- .multilpa_resample_groups(d3, rows, "id", draw)
  f <- tryCatch(suppressWarnings(.lta_refit(o, rs, 1L, 5000L, 1e-10)), error = function(e) NULL)
  if (is.null(f)) return("refit error")
  a <- tryCatch({ al <- .lta_align_estimates(f, o); sprintf("aligned %d names, setequal %s", length(al), setequal(names(al), names(.lta_bootstrap_estimates(o)))) },
                error = function(e) paste("ALIGN ERROR:", conditionMessage(e), "| class", class(e)[1]))
  sprintf("has long person %s; n_occasions %d; %s", any(draw == 1L), f$n_occasions, a)
}, "")
print(res)
for (it in c(20, 60)) {
  r <- tryCatch(withCallingHandlers({ b <- parameter_inference(o, data = d3, method = "bootstrap", iter = it, n_starts = 1, seed = 1); sprintf("completed n_valid %d of %d", attr(b, "n_valid"), it) }, warning = function(c) invokeRestart("muffleWarning")),
     error = function(e) sprintf("ERROR [%s]: %s", paste(class(e), collapse = ","), conditionMessage(e)))
  cat("iter", it, "->", r, "\n")
}
try_boot <- function(label, f, data, it = 30) {
  r <- tryCatch(withCallingHandlers({ b <- parameter_inference(f, data = data, method = "bootstrap", iter = it, n_starts = 1, seed = 2)
     sprintf("completed, n_valid %d of %d; messages: %s", attr(b, "n_valid"), it, paste(unique(na.omit(attr(b, "messages"))), collapse = " | ")) }, warning = function(c) invokeRestart("muffleWarning")),
     error = function(e) sprintf("ERROR class [%s]: %s", paste(class(e), collapse = ","), conditionMessage(e)))
  cat(label, "->", r, "\n")
}
m <- suppressWarnings(lta(d3, vars, "id", 2, "time", measurement = "occasion", n_starts = 2, seed = 1)); try_boot("occasion measurement, one longer person", m, d3)
hc <- suppressWarnings(lta(d3, vars, "id", 2, "time", n_starts = 2, seed = 1)); try_boot("homogeneous engine, one longer person", hc, d3)
rare <- sim(60, 3, 8); rare$f <- factor(ifelse(rare$id == 1, "rare", ifelse(rare$id %% 2 == 0, "a", "b")))
c1 <- suppressWarnings(lta(rare, vars, "id", 2, "time", transition_covariates = "f", n_starts = 2, seed = 1)); try_boot("factor transition covariate with a rare level", c1, rare)
c2 <- suppressWarnings(lta(rare, vars, "id", 2, "time", initial_covariates = "f", n_starts = 2, seed = 1)); try_boot("factor initial covariate with a rare level", c2, rare)
