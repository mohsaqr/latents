suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(3)
quiet <- function(e) withCallingHandlers(e, warning=function(w) {invokeRestart("muffleWarning")})
J <- 50; nj <- sample(3:6, J, TRUE); gid <- rep(seq_len(J), nj); n <- length(gid)
gc <- rbinom(J,1,.5)[gid]; z <- rbinom(n,1,ifelse(gc==1,.8,.2))+1L
dat <- data.frame(g=sample(letters,1)|>paste0(gid), y1=rnorm(n,c(0,2.5)[z],c(1,.7)[z]), y2=rnorm(n,c(1,-1.5)[z],c(.8,1.2)[z]),
  y3=rnorm(n,c(0,1)[z],1), cat=vapply(z,\(k) sample(1:3,1,prob=if(k==1) c(.7,.2,.1) else c(.1,.3,.6)),1L),
  ord=vapply(z,\(k) sample(1:4,1,prob=if(k==1) c(.5,.3,.15,.05) else c(.05,.15,.3,.5)),1L),
  cnt=rnbinom(n,size=2,mu=c(2,6)[z]), cov=rnorm(n)+z)
dat$w <- c(1,2,3)[(gid%%3)+1]
# row order shuffled so that group rows interleave
dat <- dat[sample(n),]; rownames(dat) <- NULL
miss <- dat; for (v in c("y1","y2","cat","ord","cnt")) miss[[v]][sample(n, round(.08*n))] <- NA
chk <- function(label, fit, nd=NULL) {
  p <- quiet(predict(fit, type="posterior"))
  cols <- grep("^probability_profile_", names(p), value=TRUE)
  post <- as.matrix(p[cols]); sp <- unname(fit$subject_posteriors)
  dimnames(post) <- NULL
  ok <- !is.na(post[,1])
  dmax <- max(abs(post[ok,,drop=FALSE]-sp[ok,,drop=FALSE]))
  cl <- predict(fit, type="class")
  dcl <- sum(cl$profile[ok]!=fit$subject_profiles[ok])
  gd <- NA
  if (!is.null(p$group_posterior)) {
    gp <- fit$group_posteriors
    gd <- max(abs(p$group_posterior[!is.na(p$group_posterior)] - gp[cbind(fit$group_index, p$group_class)][!is.na(p$group_posterior)]))
  }
  # density: single level sums to LL
  ds <- quiet(predict(fit, type="density"))
  data.frame(case=label, dmax_post=dmax, class_diffs=dcl, dmax_group=gd, excluded=sum(!ok),
             rowsum1=max(abs(rowSums(post[ok,,drop=FALSE]) + (if(isTRUE(fit$noise)) p$probability_noise[ok] else 0) -1)),
             density_sum=sum(ds$log_density), LL=fit$log_likelihood)
}
base <- list(data=miss, vars=c("y1","y2","y3"), id="g", n_profiles=2, n_group_classes=2, n_starts=3, seed=1, missing="fiml", tol=1e-9)
run <- function(...) quiet(do.call(multilpa, modifyList(base, list(...))))
out <- list(
 chk("VVI fiml 2gc", run()),
 chk("VVV fiml 2gc", run(model="VVV")),
 chk("VEV fiml", run(model="VEV")),
 chk("mixed nb fiml 2gc", run(vars=c("y1","y2","cat","ord","cnt"),categorical="cat",ordinal="ord",count="cnt",count_model="negative_binomial")),
 chk("weights", run(weights="w")),
 chk("direct", run(direct_effects="y1 ~ cov")),
 chk("ld cont", run(local_dependence="y1:y2")),
 chk("ld cat", run(vars=c("y1","cat","ord"),categorical=c("cat","ord"),local_dependence="cat:ord")),
 chk("grand centre", run(centering="grand")),
 chk("grand scale", run(scaling="grand")),
 chk("grand both", run(centering="grand",scaling="grand")),
 chk("person centre", run(centering="person")),
 chk("person both", run(centering="person",scaling="person")),
 chk("person centre grand scale", run(centering="person",scaling="grand")),
 chk("grand centre person scale", run(centering="grand",scaling="person")),
 chk("single", quiet(lpa(miss,c("y1","y2","y3"),2,missing="fiml",n_starts=3,seed=1))),
 chk("single noise", quiet(lpa(dat,c("y1","y2","y3"),2,noise=TRUE,n_starts=3,seed=1))),
 chk("single noise grand scale", quiet(lpa(dat,c("y1","y2","y3"),2,noise=TRUE,scaling="grand",n_starts=3,seed=1))),
 chk("prior", quiet(lpa(dat,c("y1","y2","y3"),2,prior=prior_control(),n_starts=3,seed=1))),
 chk("lca", quiet(lca(miss,c("cat","ord"),2,missing="fiml",n_starts=3,seed=1)))
)
print(do.call(rbind,out),digits=5)
# newdata equals training data; row permutation of newdata
f <- run(model="VVI")
pa <- predict(f, newdata=miss, type="posterior")
pb <- predict(f, type="posterior")
cat("newdata=train vs none diff", max(abs(as.matrix(pa[grep("prob",names(pa))])-as.matrix(pb[grep("prob",names(pb))]))), "\n")
perm <- sample(nrow(miss)); pc <- predict(f, newdata=miss[perm,], type="posterior")
cat("row permuted newdata diff", max(abs(as.matrix(pc[grep("prob",names(pc))])-as.matrix(pa[grep("prob",names(pa))])[perm,])), "\n")
# group with zero weight: predict posterior vs fit
miss$w0 <- ifelse(miss$g==miss$g[1],0,1)
fz <- run(weights="w0"); print(chk("zero weight group", fz))
# Person-scaled fit with zero-weight group degenerate: does predict(fit) work?
m2 <- miss; m2$y1[m2$g==m2$g[1]] <- 5
fz2 <- tryCatch(run(data=m2, weights="w0", scaling="person"), error=function(e) e)
if (!inherits(fz2,"error")) { r <- tryCatch(predict(fz2,type="posterior"), error=function(e) e); cat("zero-weight degenerate person scale predict:", if (inherits(r,"error")) paste("ERROR", class(r)[1], conditionMessage(r)) else "ok","\n") } else cat("fit error", conditionMessage(fz2),"\n")
# zero weight group with grand: predict offsets are the fit's
fz3 <- run(weights="w0", centering="grand", scaling="grand"); print(chk("zero weight grand",fz3))
fz4 <- run(weights="w0", centering="person", scaling="person"); print(chk("zero weight person",fz4))
