suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Hand computation, base R only, of every criterion on two-level and single-level fits
ng <- 40; nper <- 6
sim <- function() {
  g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each = nper)
  p <- ifelse(g[id] == 1, .75, .25)
  k <- 1 + rbinom(length(id), 1, p)
  data.frame(id = id, y1 = rnorm(length(id), c(-1.2, 1.2)[k], c(1, .7)[k]),
             y2 = rnorm(length(id), c(0, 1.5)[k], c(.8, 1.2)[k]),
             z = rnorm(length(id)), gz = rnorm(ng)[id])
}
d <- sim()
hand <- function(fit, d, vars = c("y1", "y2"), wt = NULL) {
  K <- fit$n_profiles; G <- fit$n_group_classes %||% 1L
  X <- as.matrix(d[vars]); idx <- match(d$id, unique(d$id))
  fk <- sapply(seq_len(K), function(k) {
    Reduce(`*`, lapply(seq_along(vars), function(j) dnorm(X[, j], fit$means[k, j], sqrt(fit$variances[k, j]))))
  })
  lj <- sapply(seq_len(G), function(h) {
    mix <- as.vector(fk %*% fit$profile_probabilities[h, ])
    as.vector(rowsum(log(mix), idx, reorder = TRUE)) + log(fit$group_probabilities[h])
  })
  lj <- matrix(lj, ncol = G)
  m <- apply(lj, 1, max); LL_j <- m + log(rowSums(exp(lj - m)))
  w <- wt %||% rep(1, length(LL_j))
  LL <- sum(w * LL_j)
  om <- exp(lj - LL_j)
  # complete-data entropy
  EN_g <- -sum(w * rowSums(ifelse(om > 0, om * log(om), 0)))
  EN_r <- 0
  for (h in seq_len(G)) { # hand check; loop is fine in a probe
    den <- as.vector(fk %*% fit$profile_probabilities[h, ])
    ph <- sweep(fk, 2, fit$profile_probabilities[h, ], "*") / den
    EN_r <- EN_r + sum(w[idx] * om[idx, h] * (-rowSums(ifelse(ph > 0, ph * log(ph), 0))))
  }
  list(LL = LL, EN = EN_g + EN_r, om = om)
}
ic <- function(LL, q, ng, nobs, EN) {
  c(aic = -2*LL + 2*q, kic = -2*LL + 3*(q+1),
    bic_groups = -2*LL + q*log(ng), bic_individual = -2*LL + q*log(nobs),
    sabic_groups = -2*LL + q*log((ng+2)/24), sabic_individual = -2*LL + q*log((nobs+2)/24),
    caic_groups = -2*LL + q*(log(ng)+1), caic_individual = -2*LL + q*(log(nobs)+1),
    awe_groups = -2*LL + 2*EN + 2*q*(1.5+log(ng)), awe_individual = -2*LL + 2*EN + 2*q*(1.5+log(nobs)),
    icl_groups = -2*LL + q*log(ng) + 2*EN, icl_individual = -2*LL + q*log(nobs) + 2*EN,
    clc_groups = -2*LL + 2*EN, clc_individual = -2*LL + 2*EN)
}
f22 <- multilpa(d, c("y1","y2"), "id", n_profiles = 2, n_group_classes = 2, model = "VVI", n_starts = 5, seed = 1)
h <- hand(f22, d)
q <- 2*2 + 2*2 + (2-1)*2 + 1
cat("VVI 2x2: LL pkg", f22$log_likelihood, "hand", h$LL, " q pkg", f22$n_parameters, "hand", q, "\n")
wide <- get_results(f22, "information_criteria", format = "wide")
hic <- ic(h$LL, q, ng, nrow(d), h$EN)
print(round(unlist(wide[names(hic)]) - hic, 8))
cat("fit$bic", f22$bic, "bic_groups", wide$bic_groups, "bic_individual", wide$bic_individual, "fit$bic_individual", f22$bic_individual, "\n")
# relative entropy at each level
e <- get_results(f22, "entropy")
print(e)
pr <- f22$subject_posteriors
cat("hand rel entropy indiv", 1 - (-sum(ifelse(pr>0,pr*log(pr),0)))/(nrow(pr)*log(2)),
    " groups", 1 - (-sum(ifelse(h$om>0,h$om*log(h$om),0)))/(ng*log(2)), "\n")
# Check against mover covariates: complete entropy on a covariate fit
fc <- multilpa(d, c("y1","y2"), "id", n_profiles = 2, n_group_classes = 2, model = "VVI",
               profile_covariates = "z", n_starts = 3, seed = 1)
cat("class:", class(fc), "\n")
wc <- get_results(fc, "information_criteria", format = "wide")
print(round(unlist(wc[c("aic","bic_groups","icl_groups","icl_individual","clc_groups","awe_groups")]),3))
cat("classification_entropy_rows stored:", !is.null(fc$classification_entropy_rows), "\n")
# hand EN for covariate fit
if (!is.null(fc$classification_entropy_rows)) {
  X <- as.matrix(d[c("y1","y2")]); idx <- match(d$id, unique(d$id))
  fk <- sapply(1:2, function(k) dnorm(X[,1], fc$means[k,1], sqrt(fc$variances[k,1])) * dnorm(X[,2], fc$means[k,2], sqrt(fc$variances[k,2])))
  str(fc$profile_design[[1]][1:2,]); print(fc$profile_coefficients); print(fc$group_coefficients)
  pri <- lapply(fc$profile_design, function(D) .multilpa_softmax(D, fc$profile_coefficients))
  lj <- sapply(1:2, function(h) as.vector(rowsum(log(rowSums(fk * pri[[h]])), idx, reorder = TRUE)))
  gp <- .multilpa_softmax(fc$group_design, fc$group_coefficients)
  lj2 <- lj + log(gp); m <- apply(lj2,1,max); LLj <- m + log(rowSums(exp(lj2-m))); om <- exp(lj2-LLj)
  cat("cov LL", fc$log_likelihood, "hand", sum(LLj), "\n")
  ENg <- -sum(om*log(om)); ENr <- sum(sapply(1:2, function(h) { ph <- fk*pri[[h]]/rowSums(fk*pri[[h]]); sum(om[idx,h] * -rowSums(ifelse(ph>0,ph*log(ph),0))) }))
  cat("cov EN pkg", .multilpa_complete_entropy(fc), " hand", ENg+ENr, "\n")
}
# weighted: integer weights == duplicated groups
gw <- sample(1:3, ng, TRUE); d$w <- gw[d$id]
fw <- multilpa(d, c("y1","y2"), "id", n_profiles = 2, n_group_classes = 2, model = "VVI", weights = "w", n_starts = 5, seed = 1)
dd <- do.call(rbind, lapply(seq_len(ng), function(j) { r <- d[d$id == j, ]; do.call(rbind, lapply(seq_len(gw[j]), function(c) { r$id <- j*100 + c; r })) }))
fd <- multilpa(dd, c("y1","y2"), "id", n_profiles = 2, n_group_classes = 2, model = "VVI", n_starts = 5, seed = 1)
cat("weights norm?", head(unname(fw$sampling_weights)), " LL weighted", fw$log_likelihood, " LL duplicated", fd$log_likelihood, "\n")
ww <- get_results(fw, "information_criteria", format = "wide"); wd <- get_results(fd, "information_criteria", format = "wide")
print(rbind(weighted = unlist(ww), duplicated = unlist(wd)))
print(rbind(get_results(fw,"entropy"), get_results(fd,"entropy")))
print(get_results(fw,"classification")); print(get_results(fd,"classification"))
cat("n sampling sum", sum(fw$sampling_weights), "\n")
# Noise/one-class/single-level
s1 <- lpa(d, c("y1","y2"), n_profiles = 2, model = "VVI", n_starts = 3, seed = 1)
ws <- get_results(s1, "information_criteria", format = "wide")
print(ws)
print(get_results(s1, "entropy"))
# LCA and LTA classes
print(class(s1))
