suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)
w <- function(expr) { warns <- character(); v <- withCallingHandlers(tryCatch(expr,error=function(e) { warns <<- c(warns,paste0("ERROR:",class(e)[1])); NULL }),warning=function(w){ warns <<- c(warns,class(w)[1]); invokeRestart("muffleWarning") }); list(v=v,w=warns) }
n <- 400; tr <- 1L+(runif(n)<.5)
d <- data.frame(a=rnorm(n,c(-2.5,2.5)[tr]),b=rnorm(n,c(-2.5,2.5)[tr]))
f <- lpa(d,c("a","b"),n_profiles=2,n_starts=3,seed=1)
d$x_strong <- ifelse(f$subject_profiles==1,1,-1)*runif(n,.2,2)
d$x_mod <- ifelse(f$subject_profiles==1,1,-1)*.6+rnorm(n,0,.15)
ro <- w(r3step(f,d,"x_mod")); ob <- w(r3step(f,d,"x_mod",vcov_type="observed"))
show <- function(a,l){ cat(l,"conditions:",a$w,"\n"); if(!is.null(a$v)) print(as.data.frame(a$v)[,c("term","estimate","standard_error","p_value")]) }
show(ro,"robust"); show(ob,"observed")
invisible(lapply(c(.1,.2,.3,.45),function(s){ set.seed(round(s*100)); d$x <- ifelse(f$subject_profiles==1,1,-1)*.6+rnorm(n,0,s); show(w(r3step(f,d,"x")),sprintf("noise sd %.2f robust",s)); show(w(r3step(f,d,"x",vcov_type="observed")),"  observed"); os <- w(multilpa(transform(d,g=seq_len(n)),c("a","b"),"g",n_profiles=2,n_group_classes=1,profile_covariates="x",n_starts=3,seed=1)); cat("  one-step fit:",os$w,"\n") }))
# multiplicity: BH invariance and the other methods, 3 classes binary
set.seed(1); N <- 600; tr3 <- sample(1:3,N,TRUE); d3 <- data.frame(a=rnorm(N,c(-2,0,2)[tr3]),b=rnorm(N,c(2,-2,2)[tr3]))
d3$yb <- factor(ifelse(runif(N)<c(.5,.55,.6)[tr3],"y","n")); d3$yn <- as.numeric(d3$yb=="y")
f3 <- lpa(d3,c("a","b"),n_profiles=3,n_starts=3,seed=1)
invisible(lapply(c("BH","holm","hochberg","hommel","bonferroni","BY"),function(m){
  tc <- as.data.frame(three_step(f3,d3,"yb",contrast="pairs",adjust=m)); tc <- subset(tc,category=="y"); tn <- as.data.frame(three_step(f3,d3,"yn",contrast="pairs",adjust=m))
  cat(sprintf("%-11s factor %s | numeric %s\n",m,paste(format(tc$p_adjusted,digits=4),collapse=" "),paste(format(tn$p_adjusted,digits=4),collapse=" ")))
}))
