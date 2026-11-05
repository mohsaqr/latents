suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
# Occasion-specific measurement: separation falls over occasions; a single pooled D is applied to all.
n <- 500; TT <- 4; z <- matrix(rnorm(n*TT),n,TT)
s <- matrix(0L,n,TT); s[,1] <- 1L+(runif(n)<.5)
invisible(lapply(2:TT,function(t){ mv <- runif(n)<plogis(-1.2+.8*z[,t]); s[,t] <<- ifelse(mv,3L-s[,t-1],s[,t-1]) }))
sep <- c(2.2,1.6,1,.6)
long <- data.frame(id=rep(seq_len(n),TT),time=rep(seq_len(TT),each=n),z=as.vector(z))
sv <- as.vector(s); st <- rep(sep,each=n)
long$a <- rnorm(n*TT,ifelse(sv==1,-st,st)); long$b <- rnorm(n*TT,ifelse(sv==1,-st,st))
args <- formals(lta); print(names(args))
fit <- lta(long,c("a","b"),id="id",time="time",n_profiles=2,measurement="occasion",n_starts=5,seed=1)
print(class(fit))
r <- r3step(fit,long,transition_covariates="z")
print(as.data.frame(r)[,c("block","from","to","term","estimate","standard_error")]); print(attr(r,"classification_errors"))
post <- .lta_occasion_posteriors(fit); W <- max.col(post,"first")
invisible(lapply(1:TT,function(t){ rows <- long$time==t; P <- post[rows,]; w <- W[rows]
  D <- t(vapply(1:2,function(k) vapply(1:2,function(j) sum(P[w==j,k])/sum(P[,k]),0),c(0,0))); cat("t",t,"D diag",format(diag(D),digits=3),"\n") }))
# Truth-side: fit with occasion-specific D by hand (independent ML)
Dl <- lapply(1:TT,function(t){ rows <- long$time==t; P <- post[rows,]; w <- W[rows]; t(vapply(1:2,function(k) vapply(1:2,function(j) sum(P[w==j,k])/sum(P[,k]),0),c(0,0))) })
M <- matrix(W,n,TT); Z <- z
nll <- function(par,Dfun){ -sum(vapply(seq_len(n),function(i){ p1 <- plogis(par[1]); a <- c(p1,1-p1)*Dfun(1)[,M[i,1]]
   for (t in 2:TT){ a12 <- plogis(par[2]+par[3]*Z[i,t]); a21 <- plogis(par[4]+par[5]*Z[i,t]); A <- rbind(c(1-a12,a12),c(a21,1-a21)); a <- as.vector(a%*%A)*Dfun(t)[,M[i,t]] }
   log(sum(a)) },0)) }
pooled <- attr(r,"classification_errors")
o1 <- optim(r$estimate,nll,Dfun=function(t) pooled,method="BFGS",control=list(reltol=1e-9))
o2 <- optim(r$estimate,nll,Dfun=function(t) Dl[[t]],method="BFGS",control=list(reltol=1e-9))
cat("pooled D (matches package):",format(o1$par,digits=4),"\noccasion D:                ",format(o2$par,digits=4),"\ntruth labels may swap; truth int -1.2 slope .8\n")
print(as.data.frame(get_results(fit,"transitions")))
