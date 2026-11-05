suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(2)
quiet <- function(e) withCallingHandlers(e, warning=function(w) {invokeRestart("muffleWarning")})
J <- 60; nj <- sample(3:6, J, TRUE); gid <- rep(seq_len(J), nj); n <- length(gid)
gc <- rbinom(J,1,.5)[gid]; z <- rbinom(n,1,ifelse(gc==1,.8,.2))+1L
dat <- data.frame(g=gid, y1=rnorm(n,c(0,2.5)[z],c(1,.7)[z]), y2=rnorm(n,c(1,-1.5)[z],c(.8,1.2)[z]),
  y3=rnorm(n,c(0,1)[z],1), cat=vapply(z,\(k) sample(1:3,1,prob=if(k==1) c(.7,.2,.1) else c(.1,.3,.6)),1L),
  ord=vapply(z,\(k) sample(1:4,1,prob=if(k==1) c(.5,.3,.15,.05) else c(.05,.15,.3,.5)),1L),
  cnt=rnbinom(n,size=2,mu=c(2,6)[z]), cov=rnorm(n))
wg <- sample(1:3, J, TRUE, prob=c(.5,.3,.2))
dat$w <- wg[gid]
miss <- dat; for (v in c("y1","y2","cat","ord","cnt")) miss[[v]][sample(n, round(.08*n))] <- NA
# Duplicate the groups: copy g, renamed
expand <- function(d) do.call(rbind, lapply(seq_len(J), function(j) {
  part <- d[d$g==j,]; do.call(rbind, lapply(seq_len(wg[j]), function(c) { part$g <- paste(j,c,sep="_"); part }))
}))
dup <- expand(miss); dup$w <- NULL
cmp <- function(label, args, single=FALSE) {
  a_w <- c(list(data=if (single) subset(miss, !duplicated(g)) else miss, weights="w"), args)
  a_d <- c(list(data=if (single) subset(expand(miss), !duplicated(g)) else dup), args)
  if (!single) { a_w$id <- "g"; a_d$id <- "g" } else { a_w$id <- NULL; a_d$id <- NULL }
  fw <- quiet(do.call(multilpa, a_w)); fd <- quiet(do.call(multilpa, a_d))
  # Align by first mean of y1 or location
  key <- function(f) if (!is.null(f$means) && ncol(f$means)>0) f$means[,1] else f$response_probabilities[[1]][,1]
  o1 <- order(key(fw)); o2 <- order(key(fd))
  sumw <- if (single) sum(subset(miss,!duplicated(g))$w) else sum(wg)
  nunits <- if (single) nrow(subset(miss,!duplicated(g))) else J
  data.frame(case=label,
    dLL_scaled = fw$log_likelihood - fd$log_likelihood*nunits/sumw,
    dMean = if (!is.null(fw$means) && length(fw$means)) max(abs(fw$means[o1,,drop=FALSE]-fd$means[o2,,drop=FALSE])) else NA,
    dVar = if (!is.null(fw$variances) && length(fw$variances)) max(abs(fw$variances[o1,,drop=FALSE]-fd$variances[o2,,drop=FALSE])) else NA,
    dPi = max(abs(colSums(fw$group_probabilities*fw$profile_probabilities)[o1]-colSums(fd$group_probabilities*fd$profile_probabilities)[o2])),
    np=fw$n_parameters, npd=fd$n_parameters, conv=fw$converged&&fd$converged, bicw=round(fw$bic,3))
}
base <- list(vars=c("y1","y2","y3"), n_profiles=2, n_group_classes=1, n_starts=6, seed=1, tol=1e-11, max_iter=20000)
out <- list(
 cmp("VVI complete", modifyList(base,list(n_group_classes=1))|>c(list(missing="fiml"))),
 cmp("EEI fiml", c(base,list(missing="fiml",model="EEI"))),
 cmp("EVI fiml", c(base,list(missing="fiml",model="EVI"))),
 cmp("VEV fiml", c(base,list(missing="fiml",model="VEV"))),
 cmp("EVE fiml", c(base,list(missing="fiml",model="EVE"))),
 cmp("VVV 2 grp classes", c(modifyList(base,list(n_group_classes=2)),list(missing="fiml",model="VVV"))),
 cmp("mixed cat+ord+nb 2gc", c(modifyList(base,list(n_group_classes=2,vars=c("y1","y2","cat","ord","cnt"))),list(missing="fiml",categorical="cat",ordinal="ord",count="cnt",count_model="negative_binomial"))),
 cmp("poisson", c(modifyList(base,list(vars=c("y1","cat","cnt"))),list(missing="fiml",categorical="cat",count="cnt"))),
 cmp("direct effects", c(base,list(missing="fiml",direct_effects="y1 ~ cov"))),
 cmp("local dep cont", c(base,list(missing="fiml",local_dependence="y1:y2"))),
 cmp("local dep cat", c(modifyList(base,list(vars=c("y1","cat","ord"))),list(missing="fiml",categorical=c("cat","ord"),local_dependence="cat:ord")))
)
print(do.call(rbind,out), digits=4)
# Scaling and centring under weights: is the reported solution the same as the duplicated one?
sc <- function(centering, scaling) {
  fw <- quiet(multilpa(miss, c("y1","y2"), "g", n_profiles=2, n_group_classes=1, weights="w", missing="fiml", centering=centering, scaling=scaling, n_starts=4, seed=1, tol=1e-10))
  fd <- quiet(multilpa(dup, c("y1","y2"), "g", n_profiles=2, n_group_classes=1, missing="fiml", centering=centering, scaling=scaling, n_starts=4, seed=1, tol=1e-10))
  o1 <- order(fw$means[,1]); o2 <- order(fd$means[,1])
  data.frame(centering, scaling, dMean=max(abs(fw$means[o1,]-fd$means[o2,])), dVar=max(abs(fw$variances[o1,]-fd$variances[o2,])),
    off_w=if (!is.null(fw$centering_offsets)&&nrow(fw$centering_offsets)==1) fw$centering_offsets[1,1] else NA,
    off_d=if (!is.null(fd$centering_offsets)&&nrow(fd$centering_offsets)==1) fd$centering_offsets[1,1] else NA,
    sc_w=if (!is.null(fw$scaling_factors)&&nrow(fw$scaling_factors)==1) fw$scaling_factors[1,1] else NA,
    sc_d=if (!is.null(fd$scaling_factors)&&nrow(fd$scaling_factors)==1) fd$scaling_factors[1,1] else NA,
    wmean_of_means=sum(colSums(fw$group_probabilities*fw$profile_probabilities)*fw$means[,1]))
}
print(do.call(rbind, list(sc("grand","none"),sc("none","grand"),sc("grand","grand"),sc("person","none"),sc("person","person"))), digits=5)
# min_variance under weights: rows vs duplicated (default bound)
fw <- quiet(multilpa(miss, c("y1","y2"), "g", n_profiles=2, n_group_classes=1, weights="w", missing="fiml", n_starts=2, seed=1))
fd <- quiet(multilpa(dup, c("y1","y2"), "g", n_profiles=2, n_group_classes=1, missing="fiml", n_starts=2, seed=1))
cat("min_variance w", fw$min_variance, " dup", fd$min_variance, "\n")
