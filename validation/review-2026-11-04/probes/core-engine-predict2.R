suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(3)
quiet <- function(e) withCallingHandlers(e, warning=function(w) {invokeRestart("muffleWarning")})
J <- 50; nj <- sample(3:6, J, TRUE); gid <- rep(seq_len(J), nj); n <- length(gid)
gc <- rbinom(J,1,.5)[gid]; z <- rbinom(n,1,ifelse(gc==1,.8,.2))+1L
dat <- data.frame(g=paste0("s",gid), y1=rnorm(n,c(0,2.5)[z],c(1,.7)[z]), y2=rnorm(n,c(1,-1.5)[z],c(.8,1.2)[z]))
dat$w0 <- ifelse(dat$g=="s1",0,1)
post_diff <- function(fit) { p <- predict(fit,type="posterior"); cols<-grep("^probability_profile_",names(p),value=TRUE); m<-as.matrix(p[cols]); dimnames(m)<-NULL; max(abs(m-unname(fit$subject_posteriors)),na.rm=TRUE) }
f<-function(...) quiet(multilpa(dat,c("y1","y2"),"g",2,2,n_starts=2,seed=1,weights="w0",...))
cat("zero weight, none", post_diff(f()),"\n")
cat("zero weight, grand both", post_diff(f(centering="grand",scaling="grand")),"\n")
cat("zero weight, person both", post_diff(f(centering="person",scaling="person")),"\n")
d2 <- dat; d2$y1[d2$g=="s1"] <- 5
ff <- quiet(multilpa(d2,c("y1","y2"),"g",2,2,n_starts=2,seed=1,weights="w0",scaling="person"))
cat("fit ok; factors for s1 y1:", ff$scaling_factors[match("s1",ff$group_ids),"y1"],"\n")
r <- tryCatch(predict(ff,type="posterior"),error=function(e) e)
cat("predict(fit) on it:", if (inherits(r,"error")) paste("ERROR:",class(r)[1],"-",conditionMessage(r)) else paste("ok, diff", post_diff(ff)),"\n")
r2 <- tryCatch(predict(ff,newdata=subset(d2,g!="s1"),type="posterior"),error=function(e) e)
cat("predict newdata without s1:", if (inherits(r2,"error")) paste("ERROR:",class(r2)[1]) else "ok","\n")
# Is a single-row new unit under person centring given a class silently?
fp <- quiet(multilpa(dat,c("y1","y2"),"g",2,2,n_starts=2,seed=1,centering="person"))
print(predict(fp,newdata=data.frame(g="new",y1=99,y2=-99),type="posterior"))
