suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(4)
quiet <- function(e) withCallingHandlers(e, warning=function(w) {invokeRestart("muffleWarning")})
n <- 200; z <- rbinom(n,1,.45)+1L
d0 <- data.frame(row=seq_len(n), y1=rnorm(n,c(0,2)[z],1), y2=NA_real_, y3=rnorm(n,c(1,-1)[z],c(.7,1.1)[z]), cv=rnorm(n)+z/2)
d0$y2 <- .6*d0$y1 + rnorm(n,c(.5,0)[z],.8)
d0$y1 <- d0$y1 + .5*d0$cv
miss <- d0; for (v in c("y1","y2","y3")) miss[[v]][sample(n,25)] <- NA
miss <- subset(miss, !is.na(y1)|!is.na(y2)|!is.na(y3))
Y <- as.matrix(miss[c("y1","y2","y3")])
ld <- function(row,mu,S){o<-which(!is.na(row)); if(!length(o)) return(0); R<-chol(S[o,o,drop=FALSE]); e<-backsolve(R,row[o]-mu[o],transpose=TRUE); -.5*(length(o)*log(2*pi)+2*sum(log(diag(R)))+sum(e^2))}
ll <- function(p,M,Sl,Y){sum(apply(Y,1,function(r){v<-vapply(1:2,\(k)log(p[k])+ld(r,M[k,],Sl[[k]]),0);m<-max(v);m+log(sum(exp(v-m)))}))}
## (a1) continuous LD + FIML, VVI: stationarity with optim
f <- quiet(lpa(miss,c("y1","y2","y3"),2,local_dependence="y1:y2",missing="fiml",n_starts=5,seed=1,tol=1e-11,max_iter=20000))
Sl <- lapply(1:2,\(k) matrix(f$covariances[,,k],3,3))
cat("LD fiml hand", ll(drop(f$profile_probabilities),unname(f$means),Sl,Y)-f$log_likelihood, " np",f$n_parameters,"\n")
print(round(Sl[[1]],4))
theta0 <- c(qlogis(f$profile_probabilities[1,1]),c(f$means),log(f$variances[1,]),log(f$variances[2,]),Sl[[1]][1,2],Sl[[2]][1,2])
obj <- function(t){p<-plogis(t[1]);M<-matrix(t[2:7],2,3);S<-lapply(1:2,function(k){s<-diag(exp(t[7+(k-1)*3+1:3]));s[1,2]<-s[2,1]<-t[13+k];s});tryCatch(-ll(c(p,1-p),M,S,Y),error=function(e) 1e10)}
op <- optim(theta0,obj,method="BFGS",control=list(maxit=200,reltol=1e-12))
cat("LD fiml optim gain",-op$value-(-obj(theta0)),"  sanity: from perturbed start recovers?", {o2<-optim(theta0+rnorm(15,0,.01),obj,method="BFGS",control=list(maxit=500,reltol=1e-13));-o2$value-(-obj(theta0))},"\n")
## (a2) direct effects + FIML, EEI and VVI
for (m in c("VVI","EEI")) {
  fd <- quiet(multilpa(miss,c("y1","y2","y3"),id=NULL,n_profiles=2,variance_model=if(m=="VVI")"varying" else "equal",direct_effects="y1 ~ cv",missing="fiml",n_starts=5,seed=1,tol=1e-11,max_iter=20000))
  ll_d <- function(t) { p<-plogis(t[1]);M<-matrix(t[2:7],2,3);b<-t[8]
    V <- if(m=="VVI") lapply(1:2,\(k) diag(exp(t[8+(k-1)*3+1:3]))) else rep(list(diag(exp(t[9:11]))),2)
    Ys<-Y; Ys[,1]<-Ys[,1]-b*miss$cv; ll(c(p,1-p),M,V,Ys)}
  t0 <- c(qlogis(fd$profile_probabilities[1,1]),c(fd$means),fd$direct_effects, if(m=="VVI") c(log(fd$variances[1,]),log(fd$variances[2,])) else log(fd$variances[1,]))
  cat("direct",m," hand-fit",ll_d(t0)-fd$log_likelihood," np",fd$n_parameters,"len",length(t0)," effect",fd$direct_effects,"\n")
  o3<-optim(t0,\(t)-ll_d(t),method="BFGS",control=list(maxit=300,reltol=1e-13))
  cat("  optim gain",-o3$value-ll_d(t0),"\n")
}
## (b) categorical LD half observed: hand likelihood
m <- 300; zz <- rbinom(m,1,.5)
cd <- data.frame(a=rbinom(m,1,.2+.6*zz),b=rbinom(m,1,.3+.4*zz),c=rbinom(m,1,.2+.5*zz),e=rbinom(m,1,.7-.4*zz),f=rbinom(m,1,.6-.3*zz))
cd$b[cd$a==1 & runif(m)<.5] <- 1
cd$a[sample(m,40)] <- NA; cd$b[sample(m,40)] <- NA; cd$e[sample(m,20)] <- NA
fl <- quiet(lca(cd,names(cd),2,local_dependence="a:b",missing="fiml",n_starts=5,seed=1,tol=1e-10))
rp <- fl$response_probabilities; print(rp[["a:b"]])
hand_cat <- function(rp,pi_) { sum(vapply(seq_len(m),function(i){
    v <- vapply(1:2,function(k){ l<-log(pi_[k])
      a<-cd$a[i];b<-cd$b[i]; J<-rp[["a:b"]][k,]; nm<-colnames(rp[["a:b"]])
      fa<-sub(":.*","",nm); fb<-sub(".*:","",nm)
      keep<-(is.na(a)|fa==as.character(a))&(is.na(b)|fb==as.character(b))
      l<-l+log(sum(J[keep]))
      for(w in c("c","e","f")) if(!is.na(cd[[w]][i])) l<-l+log(rp[[w]][k,as.character(cd[[w]][i])])
      l},0); mm<-max(v); mm+log(sum(exp(v-mm)))},0))}
cat("cat LD half-observed hand",hand_cat(rp,drop(fl$profile_probabilities))-fl$log_likelihood," np",fl$n_parameters," hand 1+2*(3+1+1+1)=13\n")
## (c) evaluate-only round trip: max_iter=0 with starting_values
rt <- function(label, fit, ...) { st<-starting_values(fit); g<-quiet(do.call(multilpa,c(list(start=st,max_iter=0L),list(...))))
  cat(sprintf("%-20s LL diff %.2e post diff %.2e\n",label,g$log_likelihood-fit$log_likelihood,max(abs(unname(g$subject_posteriors)-unname(fit$subject_posteriors))))) }
J<-40;gid<-rep(seq_len(J),each=5);dd<-data.frame(g=gid,y1=rnorm(200)+rep(rnorm(J),each=5),y2=rnorm(200),k=sample(1:3,200,TRUE))
fc<-quiet(multilpa(dd,c("y1","y2","k"),"g",2,2,categorical="k",n_starts=3,seed=1)); rt("two-level mixed",fc,data=dd,vars=c("y1","y2","k"),id="g",n_profiles=2,n_group_classes=2,categorical="k")
for (sc in list(c("grand","grand"),c("person","person"),c("person","none"))) {
  fc<-quiet(multilpa(dd,c("y1","y2"),"g",2,1,centering=sc[1],scaling=sc[2],n_starts=3,seed=1)); rt(paste(sc,collapse="/"),fc,data=dd,vars=c("y1","y2"),id="g",n_profiles=2,n_group_classes=1,centering=sc[1],scaling=sc[2]) }
fc<-quiet(multilpa(dd,c("y1","y2"),"g",2,1,model="EVE",n_starts=3,seed=1)); rt("EVE",fc,data=dd,vars=c("y1","y2"),id="g",n_profiles=2,n_group_classes=1,model="EVE")
## (e) boundary flags
nb <- data.frame(ct=c(rpois(150,2),rpois(150,9)),y=rnorm(300,rep(c(0,3),each=150)))
fn <- withCallingHandlers(lpa(nb,c("y","ct"),2,count="ct",count_model="negative_binomial",n_starts=3,seed=1),warning=function(w){cat("WARN:",class(w)[1],"\n");invokeRestart("muffleWarning")})
cat("NB dispersion",fn$count_dispersion,"  fit$boundary",fn$boundary,"\n")
cat("model table boundary:", get_results(fn,"model")$boundary, "\n"); r<-tryCatch(parameter_inference(fn),error=function(e) e); cat("Wald at NB floor:", if(inherits(r,"error")) paste("REFUSED",class(r)[1]) else paste("given; rows", nrow(as.data.frame(r))), "\n"); if(!inherits(r,"error")) print(head(subset(as.data.frame(r), grepl("dispersion",parameter)),3))
lc <- data.frame(a=rep(c(1,0),each=100),b=c(rep(1,100),rbinom(100,1,.5)),c=rbinom(200,1,.5),e=rbinom(200,1,.5))
fp <- withCallingHandlers(lca(lc,c("a","b","c","e"),2,n_starts=3,seed=1),warning=function(w){cat("WARN:",class(w)[1],"\n");invokeRestart("muffleWarning")})
cat("min resp prob",min(unlist(fp$response_probabilities)),"  fit$boundary",fp$boundary,"\n"); cat("model boundary",get_results(fp,"model")$boundary,"\n"); r<-tryCatch(parameter_inference(fp),error=function(e) e); cat("Wald at prob floor:", if(inherits(r,"error")) paste("REFUSED",class(r)[1]) else "given", "\n")
## (f) drop_empty and the criteria
de <- dd; de$y1[1:3] <- NA; de$y2[1:3] <- NA; de$lone <- seq_len(nrow(de))
a<-quiet(lpa(de,c("y1","y2"),2,missing="fiml",n_starts=2,seed=1,drop_empty=TRUE)); b<-quiet(lpa(de,c("y1","y2"),2,missing="fiml",n_starts=2,seed=1,drop_empty=FALSE))
cat("single level: bic TRUE",a$bic," FALSE",b$bic,"; bic_individual",a$bic_individual,b$bic_individual,"; n rows 200, informative",a$n_informative,"\n")
## (h) unit rescale invariance
cc <- c(1000,.01); sc <- dd; sc$y1<-sc$y1*cc[1]; sc$y2<-sc$y2*cc[2]
vapply(c("EEI","VVI","EEE","VEV","EVI"),function(mm){f1<-quiet(lpa(dd,c("y1","y2"),2,model=mm,n_starts=3,seed=1,tol=1e-10));f2<-quiet(lpa(sc,c("y1","y2"),2,model=mm,n_starts=3,seed=1,tol=1e-10)); f1$log_likelihood-(f2$log_likelihood+200*sum(log(cc)))},0)|>print()
## (x) vars order and label mapping
J<-40;m<-200;zz<-rbinom(m,1,.5)+1L
vv<-data.frame(y1=rnorm(m,c(0,3)[zz]),y2=rnorm(m,c(5,-5)[zz],c(1,2)[zz]),cat=vapply(zz,\(k) sample(c("x","y"),1,prob=if(k==1)c(.9,.1) else c(.2,.8)),""),ord=vapply(zz,\(k) sample(1:3,1,prob=if(k==1)c(.8,.1,.1) else c(.1,.1,.8)),1L),ct=rpois(m,c(1,5)[zz]),cv=rnorm(m))
vv$y2<-vv$y2+2*vv$cv
A<-quiet(lpa(vv,c("y1","y2","cat","ord","ct"),2,categorical="cat",ordinal="ord",count="ct",direct_effects="y2 ~ cv",n_starts=3,seed=1,tol=1e-10))
B<-quiet(lpa(vv,c("ct","cat","y2","ord","y1"),2,categorical="cat",ordinal="ord",count="ct",direct_effects="y2 ~ cv",n_starts=3,seed=1,tol=1e-10))
cat("vars order LL diff",A$log_likelihood-B$log_likelihood," effect",A$direct_effects,B$direct_effects,"\n")
ob<-order(B$means[,"y1"]);oa<-order(A$means[,"y1"])
cat("means diff",max(abs(A$means[oa,c("y1","y2")]-B$means[ob,c("y1","y2")]))," count diff",max(abs(A$count_means[oa,]-B$count_means[ob,])),"\n")
