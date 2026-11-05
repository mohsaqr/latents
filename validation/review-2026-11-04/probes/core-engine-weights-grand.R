suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(8)
quiet <- function(e) withCallingHandlers(e, warning=function(w) invokeRestart("muffleWarning"))
# Single-level, oversampled minority profile: weights undo the oversampling.
n <- 400; z <- rbinom(n,1,.5)
d <- data.frame(a=rnorm(n,z*3,1), b=rnorm(n,-z*2,c(1,2)[z+1]), c=rnorm(n,z,1))
d$w <- ifelse(z==1,1L,3L)          # integer weights: class 0 rows count three times
dup <- d[rep(seq_len(n), d$w),]; dup$w <- NULL; rownames(dup) <- NULL
cmp <- function(model, centering, scaling) {
  fw <- quiet(lpa(d,c("a","b","c"),2,model=model,weights="w",centering=centering,scaling=scaling,n_starts=5,seed=1,tol=1e-10))
  fd <- quiet(lpa(dup,c("a","b","c"),2,model=model,centering=centering,scaling=scaling,n_starts=5,seed=1,tol=1e-10))
  o1 <- order(fw$means[,1]); o2 <- order(fd$means[,1])
  ss <- function(f) if (is.null(f$scaling_factors)) 1 else f$scaling_factors[1,]
  oo <- function(f) if (is.null(f$centering_offsets)) 0 else f$centering_offsets[1,]
  # back on the raw scale
  raw_w <- sweep(sweep(fw$means[o1,,drop=FALSE],2,ss(fw),"*"),2,oo(fw),"+"); raw_d <- sweep(sweep(fd$means[o2,,drop=FALSE],2,ss(fd),"*"),2,oo(fd),"+")
  data.frame(model,centering,scaling,reported_mean_diff=round(max(abs(fw$means[o1,]-fd$means[o2,])),5),
    raw_mean_diff=round(max(abs(raw_w-raw_d)),5), reported_var_diff=round(max(abs(fw$variances[o1,]-fd$variances[o2,])),5),
    share_diff=round(max(abs(fw$profile_probabilities[1,o1]-fd$profile_probabilities[1,o2])),5),
    wmean_profile_a=round(sum(fw$profile_probabilities[1,]*fw$means[,1]),4))
}
tab <- do.call(rbind,list(cmp("VVI","none","none"),cmp("VVI","grand","none"),cmp("VVI","grand","grand"),cmp("EII","none","none"),cmp("EII","grand","grand"),cmp("VII","grand","grand"),cmp("EEV","none","none"),cmp("EEV","grand","grand"),cmp("VEV","grand","grand")))
print(tab)
# Zero-weight row with a missing continuous value under missing = "error"
d2 <- d; d2$w <- 1L; d2$w[1:3] <- 0L
d3 <- d2; d3$a[1] <- NA
r <- tryCatch(quiet(lpa(d3,c("a","b","c"),2,weights="w",n_starts=2,seed=1)),error=function(e) e)
cat("zero-weight NA continuous, missing=error:", if (inherits(r,"error")) paste("REFUSED",class(r)[1]) else "fit","\n")
d4 <- d2; d4$k <- sample(1:2,n,TRUE); d4$k[1] <- NA
r <- tryCatch(quiet(lpa(d4,c("a","k"),2,categorical="k",weights="w",n_starts=2,seed=1)),error=function(e) e)
cat("zero-weight NA categorical, missing=error:", if (inherits(r,"error")) paste("REFUSED",class(r)[1]) else "fit","\n")
