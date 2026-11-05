suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# A handful of data sets: bounded 1-5 outcome, K = 3, moderate separation. Does the
# multi-start selection move ML away from the start Latent GOLD uses (proportional)?
em_from <- function(pieces, y, m0, v0) {
  D <- .multilpa_error_matrix(pieces); pr <- colMeans(pieces$posteriors)
  wt <- t(D[,pieces$modal]) * matrix(pr,length(y),length(pr),byrow=TRUE)
  m <- m0; v <- v0; old <- -Inf
  for (it in 1:20000) { # EM reference; each iterate needs the last
    dens <- vapply(seq_along(m),function(k) dnorm(y,m[k],sqrt(v[k])),y)
    j <- wt*dens; ll <- sum(log(rowSums(j))); r <- j/rowSums(j)
    if (abs(ll-old)<1e-11*(1+abs(ll))) break
    old <- ll; m <- as.vector(crossprod(r,y))/colSums(r)
    v <- colSums(r*outer(y,m,"-")^2)/colSums(r)
  }
  list(m=m,v=v,ll=ll)
}
one <- function(seed, sep=.8, n=400, K=3, label="") {
  set.seed(seed)
  tr <- sample.int(K,n,TRUE); mus <- list(c(-sep,-sep),c(sep,-sep),c(0,sep))
  x <- data.frame(a=rnorm(n,vapply(mus,`[`,0,1)[tr]),b=rnorm(n,vapply(mus,`[`,0,2)[tr]))
  x$y <- pmin(pmax(round(3+c(-1,0,1)[tr]+rnorm(n,0,1)),1),5)
  f <- lpa(x,c("a","b"),n_profiles=K,n_starts=5,seed=1)
  pc <- .multilpa_level_assignments(f,"individuals")
  # truth means of the rounded outcome, in fit's labels by majority
  map <- apply(table(f$subject_profiles,tr),1,which.max)
  truth <- tapply(x$y,tr,mean)[map]
  p0 <- as.vector(crossprod(pc$posteriors,x$y))/colSums(pc$posteriors)
  v0 <- colSums(pc$posteriors*outer(x$y,p0,"-")^2)/colSums(pc$posteriors)
  lg <- em_from(pc,x$y,p0,v0)
  sel <- .multilpa_step_ml_fit(pc,x$y,"varying")
  bch <- three_step(f,x,"y")
  mlt <- three_step(f,x,"y",method="ml")
  cat(sprintf("seed %d D-diag %s\n  truth %s\n  BCH   %s (se %s)\n  EM from proportional %s LL %.3f var %s\n  selected             %s LL %.3f var %s se %s\n",
    seed,paste(format(diag(.multilpa_error_matrix(pc)),digits=2),collapse=" "),
    paste(format(truth,digits=3),collapse=" "),
    paste(format(bch$estimate,digits=3),collapse=" "),paste(format(bch$standard_error,digits=2),collapse=" "),
    paste(format(lg$m,digits=3),collapse=" "),lg$ll,paste(format(lg$v,digits=3),collapse=" "),
    paste(format(sel$means,digits=3),collapse=" "),sel$log_likelihood,paste(format(sel$variances,digits=3),collapse=" "),paste(format(mlt$standard_error,digits=2),collapse=" ")))
}
invisible(lapply(1:6, one))
invisible(lapply(7:9, one, sep=1.2))
