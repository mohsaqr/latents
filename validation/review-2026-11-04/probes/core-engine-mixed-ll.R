suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
quiet <- function(e) withCallingHandlers(e, warning=function(w) {message("WARN[",class(w)[1],"]: ",conditionMessage(w)); invokeRestart("muffleWarning")})
# Two-level data with continuous, categorical, ordinal, count (Poisson, NB), missing cells.
J <- 40; nj <- sample(4:9, J, TRUE); gid <- rep(seq_len(J), nj); n <- length(gid)
gc <- rbinom(J,1,.5)[gid]; z <- rbinom(n,1,ifelse(gc==1,.8,.25))+1L
dat <- data.frame(g=paste0("s",gid),
  y1=rnorm(n,c(0,2)[z],c(1,.7)[z]), y2=rnorm(n,c(1,-1)[z],c(.8,1.2)[z]),
  cat=sample.int(3,n,TRUE)+0L, ord=NA_integer_, cnt=rpois(n,c(1,4)[z]), nb=rnbinom(n,size=2,mu=c(2,6)[z]))
dat$cat <- vapply(z,\(k) sample(1:3,1,prob=if(k==1) c(.7,.2,.1) else c(.1,.3,.6)),1L)
dat$ord <- vapply(z,\(k) sample(1:4,1,prob=if(k==1) c(.5,.3,.15,.05) else c(.05,.15,.3,.5)),1L)
dm <- dat
for (v in c("y1","y2","cat","ord","cnt","nb")) dm[[v]][sample(n, round(.1*n))] <- NA
hand <- function(fit, d, nb=FALSE, ctr=FALSE) {
  P <- fit$n_profiles; H <- fit$n_group_classes
  lev <- fit$categorical_levels$cat
  f <- matrix(0,nrow(d),P)
  for (k in seq_len(P)) {
    l <- rep(0,nrow(d))
    for (v in c("y1","y2")) { o <- !is.na(d[[v]]); l[o] <- l[o]+dnorm(d[[v]][o],fit$means[k,v],sqrt(fit$variances[k,v]),log=TRUE) }
    o <- !is.na(d$cat); l[o] <- l[o]+log(fit$response_probabilities$cat[k,match(as.character(d$cat[o]),lev)])
    # adjacent-category by hand
    a <- c(0,fit$ordinal_intercepts$ord); eta <- fit$ordinal_locations[k,"ord"]; s <- 0:3
    pr <- exp(a+eta*s); pr <- pr/sum(pr); o <- !is.na(d$ord); l[o] <- l[o]+log(pr[d$ord[o]])
    o <- !is.na(d$cnt); 
    if (nb) { r <- 1/fit$count_dispersion[k,"cnt"]; mu <- fit$count_means[k,"cnt"]; y <- d$cnt[o]
      l[o] <- l[o]+lgamma(y+r)-lgamma(r)-lgamma(y+1)+r*log(r/(r+mu))+y*log(mu/(r+mu))
    } else l[o] <- l[o]+d$cnt[o]*log(fit$count_means[k,"cnt"])-fit$count_means[k,"cnt"]-lgamma(d$cnt[o]+1)
    f[,k] <- l
  }
  gs <- split(seq_len(nrow(d)), match(d$g,unique(d$g)))
  sum(vapply(gs, function(rows) {
    lh <- vapply(seq_len(H), function(h) {
      log(fit$group_probabilities[h])+sum(vapply(rows, function(i) {
        v <- log(fit$profile_probabilities[h,])+f[i,]; m <- max(v); m+log(sum(exp(v-m))) },0))
    },0); m <- max(lh); m+log(sum(exp(lh-m)))
  },0))
}
vars <- c("y1","y2","cat","ord","cnt")
f1 <- quiet(multilpa(dm, vars, "g", n_profiles=2, n_group_classes=2, categorical="cat", ordinal="ord", count="cnt", missing="fiml", n_starts=4, seed=1, tol=1e-10))
cat("poisson LL fit", format(f1$log_likelihood,digits=12), " hand", format(hand(f1,dm),digits=12), "\n")
# n params by hand: mixing (H-1)+H(K-1)=1+2=3; means 2*2=4; var 2*2=4; cat 2*(3-1)=4; ord (4-1)+(2-1)=4; cnt 2 => 21
cat("n_par", f1$n_parameters, "hand 21\n")
f2 <- quiet(multilpa(dm, vars, "g", n_profiles=2, n_group_classes=2, categorical="cat", ordinal="ord", count="cnt", count_model="negative_binomial", missing="fiml", n_starts=4, seed=1, tol=1e-10))
cat("NB LL fit", format(f2$log_likelihood,digits=12), " hand", format(hand(f2,dm,nb=TRUE),digits=12), " n_par", f2$n_parameters, "hand 23\n")
f3 <- quiet(multilpa(dm, vars, "g", n_profiles=2, n_group_classes=2, categorical="cat", ordinal="ord", count="cnt", count_model="negative_binomial", count_dispersion="equal", missing="fiml", n_starts=4, seed=1, tol=1e-10))
cat("NB eq LL fit", format(f3$log_likelihood,digits=12), " hand", format(hand(f3,dm,nb=TRUE),digits=12), " n_par", f3$n_parameters, "hand 22\n")
# Stationarity of the ordinal+count parameters: score by central differences at the fit, hand LL
grad_check <- function(fit, nb=FALSE) {
  g <- function(set, h=1e-5) { up <- set(fit, h); dn <- set(fit,-h); (hand(up,dm,nb)-hand(dn,dm,nb))/(2*h) }
  res <- c(
   a2=g(\(f,h){f$ordinal_intercepts$ord[1] <- f$ordinal_intercepts$ord[1]+h;f}),
   a3=g(\(f,h){f$ordinal_intercepts$ord[2] <- f$ordinal_intercepts$ord[2]+h;f}),
   a4=g(\(f,h){f$ordinal_intercepts$ord[3] <- f$ordinal_intercepts$ord[3]+h;f}),
   eta1=g(\(f,h){f$ordinal_locations[1,"ord"] <- f$ordinal_locations[1,"ord"]+h;f}),
   lam1=g(\(f,h){f$count_means[1,"cnt"] <- f$count_means[1,"cnt"]+h;f}),
   lam2=g(\(f,h){f$count_means[2,"cnt"] <- f$count_means[2,"cnt"]+h;f}))
  if (nb) res <- c(res, a1=g(\(f,h){f$count_dispersion[1,"cnt"] <- f$count_dispersion[1,"cnt"]+h;f}),
     a2=g(\(f,h){f$count_dispersion[2,"cnt"] <- f$count_dispersion[2,"cnt"]+h;f}))
  res
}
print(round(grad_check(f1),5)); print(round(grad_check(f2,TRUE),5))
cat("mean sanity: fit locs", f1$ordinal_locations,"\n")
# all-ordinal LCA hand count, two ordinal 3-cat + ...; local dependence counts
ld <- data.frame(g=rep(1:60,each=4), a=sample(1:2,240,TRUE), b=sample(1:3,240,TRUE), c=sample(1:2,240,TRUE), x1=rnorm(240),x2=rnorm(240),x3=rnorm(240), w=rnorm(240))
fl <- quiet(multilpa(ld, c("a","b","c"), "g", n_profiles=2,n_group_classes=2,categorical=c("a","b","c"),local_dependence="a:b",n_starts=2,seed=1))
cat("LD cat n_par", fl$n_parameters, "hand: mix 3 + 2*(2*3-1)=10 + 2*(2-1)=2 =>15\n")
fc <- quiet(multilpa(ld, c("x1","x2","x3"), "g", n_profiles=2,n_group_classes=1,local_dependence="x1:x2",n_starts=2,seed=1))
cat("LD cont n_par VVI+cov", fc$n_parameters, "hand: mix 1 + means 6 + var 6 + cov 2 =15\n")
fc2 <- quiet(multilpa(ld, c("x1","x2","x3"), "g", n_profiles=2,n_group_classes=1,local_dependence="x1:x2",variance_model="equal",n_starts=2,seed=1))
cat("LD cont n_par EEI+cov", fc2$n_parameters, "hand: mix 1 + means 6 + var 3 + cov 1 =11\n")
fd <- quiet(multilpa(ld, c("x1","x2","x3"), "g", n_profiles=2,n_group_classes=1,direct_effects="x1 ~ w",n_starts=2,seed=1))
cat("direct n_par", fd$n_parameters, "hand: 1+6+6+1=14\n")
fn <- quiet(lpa(ld, c("x1","x2","x3"), n_profiles=2,noise=TRUE,n_starts=2,seed=1))
cat("noise n_par VVI", fn$n_parameters, "hand mclust: 1+6+6+2=15\n")
# Hand LL under noise
pi_ <- c(fn$profile_probabilities[1,], fn$noise_probability)
lln <- sum(apply(as.matrix(ld[,c("x1","x2","x3")]),1,function(r){
  v <- c(vapply(1:2,\(k) log(pi_[k])+sum(dnorm(r,fn$means[k,],sqrt(fn$variances[k,]),log=TRUE)),0), log(pi_[3])-log(fn$hypervolume))
  m <- max(v); m+log(sum(exp(v-m)))}))
cat("noise LL fit",format(fn$log_likelihood,digits=12)," hand",format(lln,digits=12)," pi sum",sum(pi_),"\n")
if (requireNamespace("mclust",quietly=TRUE)) cat("hypvol mclust", mclust::hypvol(as.matrix(ld[,c("x1","x2","x3")])), " fit", fn$hypervolume,"\n")
# Fixed: means held, hand count
ff0 <- quiet(multilpa(ld, c("x1","x2","a"), "g", n_profiles=2, n_group_classes=1, categorical="a", n_starts=2, seed=1))
st <- starting_values(ff0)
ff <- quiet(multilpa(ld, c("x1","x2","a"), "g", n_profiles=2, n_group_classes=1, categorical="a", n_starts=1, seed=1, start=st, fixed="means"))
cat("fixed means n_par",ff$n_parameters,"hand: 1+var 4 + cat 2=7; with_measurement", ff$n_parameters_with_measurement, "hand 11\n")
ff <- quiet(multilpa(ld, c("x1","x2","a"), "g", n_profiles=2, n_group_classes=1, categorical="a", n_starts=1, seed=1, start=st, fixed="measurement"))
cat("fixed meas n_par",ff$n_parameters,"hand 1; means equal", isTRUE(all.equal(unname(ff$means),unname(ff0$means))), "resp equal", isTRUE(all.equal(ff$response_probabilities$a,ff0$response_probabilities$a)),"\n")
ff <- quiet(multilpa(ld, c("x1","x2"), "g", n_profiles=2, n_group_classes=1, model="VEV", n_starts=2, seed=1))
st2 <- starting_values(ff); ff2 <- quiet(multilpa(ld, c("x1","x2"), "g", n_profiles=2, n_group_classes=1, model="VEV", n_starts=1, seed=1, start=st2, fixed="means"))
cat("VEV fixed means n_par",ff2$n_parameters," hand: VEV d=2 k=2: 2+1+2=5 + mix 1 = 6; LL gain vs free", ff$log_likelihood-ff2$log_likelihood, " boundary", ff2$boundary, "\n")
