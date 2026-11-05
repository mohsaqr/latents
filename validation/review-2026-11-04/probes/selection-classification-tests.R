suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(7)
# ---- VLMR: fully independent reference on K=1 v K=2, univariate equal variance + numDeriv
n <- 200; k <- rbinom(n,1,.4); dv <- data.frame(y = rnorm(n, ifelse(k==1, 2.2, 0)), id = seq_len(n))
grid <- enumerate_lpa(dv, "y", n_profiles = 1:2, model = "EEI", n_starts = 5, seed = 1)
lm <- get_results(grid, "lmr"); print(lm[,c("null_candidate","alternative_candidate","n_units","parameters_added","lr","vlmr_p","adjusted_lr","almr_p","reference_mean","reference_sd","note")])
f1 <- get_candidate(grid, n_profiles=1); f2 <- get_candidate(grid, n_profiles=2)
ll1 <- function(th, y=dv$y) dnorm(y, th[1], exp(th[2]), log=TRUE)
ll2 <- function(th, y=dv$y) log(plogis(th[1])*dnorm(y,th[2],exp(th[4])) + (1-plogis(th[1]))*dnorm(y,th[3],exp(th[4])))
t1 <- c(f1$means[1,1], .5*log(f1$variances[1,1]))
t2 <- c(qlogis(f2$profile_probabilities[1,1]), f2$means[1,1], f2$means[2,1], .5*log(f2$variances[1,1]))
cat("LL check", sum(ll1(t1)) - f1$log_likelihood, sum(ll2(t2)) - f2$log_likelihood, "\n")
if (requireNamespace("numDeriv", quietly=TRUE)) {
  S1 <- numDeriv::jacobian(ll1, t1); S2 <- numDeriv::jacobian(ll2, t2)
  H1 <- numDeriv::hessian(function(th) sum(ll1(th)), t1); H2 <- numDeriv::hessian(function(th) sum(ll2(th)), t2)
  B1 <- crossprod(S1); B2 <- crossprod(S2); B21 <- crossprod(S2, S1) # alt by null
  W <- rbind(cbind(-B2 %*% solve(H2), -B21 %*% solve(H1)), cbind(t(B21) %*% solve(H2), B1 %*% solve(H1)))
  ev <- Re(eigen(W)$values); print(round(ev, 6))
  lr <- 2*(sum(ll2(t2)) - sum(ll1(t1)))
  cat("reference mean hand", sum(ev), "pkg", lm$reference_mean, " sd hand", sqrt(2*sum(ev^2)), "pkg", lm$reference_sd, "\n")
  sims <- rowSums(sapply(ev[abs(ev)>1e-8], function(l) l*rchisq(2e5,1)))
  cat("MC p (2e5) VLMR", mean(sims > lr), "pkg", lm$vlmr_p, "; aLMR", mean(sims > lr/(1+1/(3*0+2*log(n)))), "pkg", lm$almr_p, "\n")
  cat("adjusted hand", lr/(1+1/(2*log(n))), "pkg", lm$adjusted_lr, "\n")
  cat("Imhof check: chisq3 at 5:", .latents_imhof_upper(5, c(1,1,1)), pchisq(5,3,lower=FALSE), "\n")
}
# ---- VLMR under weights: is B computed with weight^2, and is it refused or not?
dv$w <- sample(c(1,2,3), n, TRUE)
gw <- enumerate_lpa(dv, "y", n_profiles = 1:2, model = "EEI", n_starts = 5, seed = 1, weights = "w")
lw <- tryCatch(get_results(gw, "lmr"), error = function(e) e)
if (inherits(lw, "error")) { cat("weighted LMR:", class(lw)[1], conditionMessage(lw), "\n") } else {
  print(lw[,c("n_units","parameters_added","lr","vlmr_p","almr_p","reference_mean","reference_sd","note")])
  fw1 <- get_candidate(gw, n_profiles=1); fw2 <- get_candidate(gw, n_profiles=2)
  wn <- unname(fw1$sampling_weights); cat("sum w", sum(wn), "\n")
  th1 <- c(fw1$means[1,1], .5*log(fw1$variances[1,1])); th2 <- c(qlogis(fw2$profile_probabilities[1,1]), fw2$means[1,1], fw2$means[2,1], .5*log(fw2$variances[1,1]))
  S1 <- numDeriv::jacobian(ll1, th1)*wn; S2 <- numDeriv::jacobian(ll2, th2)*wn
  H1 <- numDeriv::hessian(function(th) sum(wn*ll1(th)), th1); H2 <- numDeriv::hessian(function(th) sum(wn*ll2(th)), th2)
  B1 <- crossprod(S1); B2 <- crossprod(S2); B21 <- crossprod(S2,S1)
  W <- rbind(cbind(-B2 %*% solve(H2), -B21 %*% solve(H1)), cbind(t(B21) %*% solve(H2), B1 %*% solve(H1)))
  ev <- Re(eigen(W)$values); cat("weighted hand mean", sum(ev), " sd", sqrt(2*sum(ev^2)), "; LL hand diff", sum(wn*ll2(th2)) - fw2$log_likelihood, "\n")
}
# ---- BLRT p-value: read replicates
dv2 <- data.frame(id = rep(1:30, each = 4)); kk <- rep(rbinom(30,1,.5), each=4); dv2$y1 <- rnorm(120, 2*kk); dv2$y2 <- rnorm(120, -2*kk)
n0 <- multilpa(dv2, c("y1","y2"), "id", n_profiles = 1, n_group_classes = 1, model = "EEI", n_starts = 2, seed = 1)
a0 <- multilpa(dv2, c("y1","y2"), "id", n_profiles = 2, n_group_classes = 1, model = "EEI", n_starts = 2, seed = 1)
b <- bootstrap_lrt(n0, a0, iter = 19, seed = 1)
rp <- get_results(b, "replicates"); print(as.data.frame(b)); print(summary(rp$statistic))
cat("hand p", (1+sum(rp$statistic[rp$valid] >= b$statistic))/(sum(rp$valid)+1), " pkg", b$p_value, " observed", 2*(a0$log_likelihood-n0$log_likelihood), "\n")
# Null near-true: replicate statistics should not exceed by order of magnitude
a1 <- multilpa(dv2, c("y1","y2"), "id", n_profiles = 3, n_group_classes = 1, model = "EEI", n_starts = 3, seed = 1)
b2 <- bootstrap_lrt(a0, a1, iter = 19, seed = 2); print(as.data.frame(b2)[,c("statistic","p_value","monte_carlo_se","n_valid","n_invalid")])
# ---- Round trip: BLRT's refit rule reproduces each spec on the original data
set.seed(8); ng <- 40; nper <- 5; g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each=nper); k <- 1+rbinom(length(id),1,ifelse(g[id]==1,.8,.2))
d <- data.frame(id=id, y1=rnorm(length(id),c(-1.5,1.5)[k]), y2=rnorm(length(id),c(0,2)[k]), z=rnorm(length(id)))
d$a <- ifelse(runif(nrow(d)) < c(.2,.8)[k], "x", "y"); d$b <- ifelse(runif(nrow(d)) < c(.7,.3)[k], "u", "v"); d$o <- factor(pmin(3,pmax(1,round(2+c(-1,1)[k]*.7+rnorm(nrow(d),0,.6)))),ordered=TRUE); d$cnt <- rnbinom(nrow(d), size=2, mu=c(1,5)[k])
specs <- list(
  VEI = list(vars=c("y1","y2"), model="VEI"), EVE = list(vars=c("y1","y2"), model="EVE"),
  pair = list(vars=c("y1","y2"), variance_model="varying", local_dependence="y1:y2"),
  direct = list(vars=c("y1","y2"), model="VVI", direct_effects="y1 ~ z"),
  joint = list(vars=c("a","b","y1"), categorical=c("a","b"), local_dependence="a:b", model="VVI"),
  ordinal = list(vars=c("o","y1"), ordinal="o", model="VVI"),
  nb = list(vars=c("cnt","y1"), count="cnt", count_model="negative_binomial", count_dispersion="equal", model="VVI"),
  grand = list(vars=c("y1","y2"), model="VVI", centering="grand"), scale = list(vars=c("y1","y2"), model="VVI", scaling="grand"))
res <- lapply(names(specs), function(nm) {
  s <- specs[[nm]]
  mk <- function(K, G) do.call(multilpa, c(list(data=d, id="id", n_profiles=K, n_group_classes=G, n_starts=2, seed=1), s))
  tryCatch({
    n <- mk(2,1); alt <- mk(2,2)
    bb <- suppressWarnings(bootstrap_lrt(n, alt, iter=3, seed=1))
    # read refit shape: replicate the BLRT refit on the data directly
    fam <- .multilpa_refit_measurement(n)
    re <- do.call(multilpa, c(list(data=d, vars=n$vars, id=n$id, n_profiles=2, n_group_classes=1), fam,
      list(n_starts=2, seed=1, min_variance=n$min_variance, centering=n$centering %||% "none", scaling=n$scaling %||% "none",
           categorical=n$categorical %||% character(), ordinal=n$ordinal %||% character(), count=n$count %||% character(),
           count_model=n$extra_data$count_model %||% "poisson", count_dispersion=n$extra_data$count_dispersion %||% "varying", min_probability=n$min_probability %||% 1e-10)))
    data.frame(spec=nm, q_fit=n$n_parameters, q_refit=if(inherits(re,"error")) NA else re$n_parameters, ll_diff=re$log_likelihood-n$log_likelihood, valid=bb$n_valid, invalid=bb$n_invalid, stat=round(bb$statistic,3), err = NA_character_)
  }, error=function(e) data.frame(spec=nm,q_fit=NA,q_refit=NA,ll_diff=NA,valid=NA,invalid=NA,stat=NA,err=substr(conditionMessage(e),1,120)))
})
print(do.call(rbind,res))
# Simulated moments: categorical/joint, NB
nj <- do.call(multilpa, c(list(data=d, id="id", n_profiles=2, n_group_classes=2, n_starts=2, seed=1), specs$joint))
sim <- do.call(rbind, lapply(1:50, function(i) .multilpa_simulate(nj)))
print(names(sim)); print(round(prop.table(table(sim[["a:b"]])),3)); 
pmix <- colSums(as.vector(nj$group_probabilities %*% nj$profile_probabilities) * nj$response_probabilities[["a:b"]]); print(round(pmix,3))
nbf <- do.call(multilpa, c(list(data=d, id="id", n_profiles=2, n_group_classes=2, n_starts=2, seed=1), specs$nb))
sim <- do.call(rbind, lapply(1:50, function(i) .multilpa_simulate(nbf)))
wk <- as.vector(nbf$group_probabilities %*% nbf$profile_probabilities); mu <- nbf$count_means[,1]; disp <- nbf$count_dispersion[1,1]
cat("NB sim mean", mean(sim$cnt), "model", sum(wk*mu), " var sim", var(sim$cnt), "model", sum(wk*(mu+disp*mu^2+mu^2))-sum(wk*mu)^2, "\n")
# ---- deltas, knees, Schwarz by hand
gg <- enumerate_lpa(dv2, c("y1","y2"), n_profiles = 1:5, model = "EEI", n_starts = 3, seed = 1)
tb <- get_results(gg, "candidates"); print(tb[,c("candidate","bic_groups","bic_individual")])
dl <- get_results(gg, "deltas"); print(dl[,c("candidate","value","delta_profile","from_profile")])
kn <- get_results(gg, "knees"); print(kn[,c("direction","knee_candidate","distance","knee_without_largest","stable")])
x <- 1:5; y <- tb$bic_groups; xs <- (x-1)/4; ys <- (y-min(y))/(max(y)-min(y)); sl <- ys[5]-ys[1]; dist <- ((ys[1]+sl*xs)-ys)/sqrt(1+sl^2); cat("hand distances", round(dist,4), "\n")
cm <- as.data.frame(compare_models(gg)); dlt <- tb$bic_groups-min(tb$bic_groups); wgt <- exp(-dlt/2)/sum(exp(-dlt/2)); cat("Schwarz diff", max(abs(cm$bic_weight-wgt)), "\n")
print(get_results(gg,"best")[,c("criterion","convention","candidate","value")])
