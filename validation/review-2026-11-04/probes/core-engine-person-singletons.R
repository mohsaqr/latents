suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(7)
caught <- function(e) { w <- character(); v <- withCallingHandlers(e, warning=function(c){w<<-c(w,class(c)[1]);invokeRestart("muffleWarning")}); attr(v,"w") <- unique(w); v }
J <- 80; sizes <- c(rep(1L,12), sample(5:9, J-12, TRUE)); gid <- rep(seq_len(J), sizes); n <- length(gid)
z <- rbinom(n,1,.5)+1L; off <- rep(rnorm(J,0,1.5),sizes)
d <- data.frame(g=gid, a=off+rnorm(n,c(-1,1)[z],.7), b=off+rnorm(n,c(1,-1)[z],.7), c=rnorm(n,c(0,.5)[z],1))
cat("rows",n,"singleton groups",sum(sizes==1),"\n")
fa <- caught(multilpa(d,c("a","b","c"),"g",n_profiles=3,n_group_classes=1,centering="person",n_starts=5,seed=1))
cat("all groups: warnings:",attr(fa,"w"),"\n"); print(round(fa$means,3)); print(round(fa$variances,5)); print(round(drop(fa$profile_probabilities),3)); cat("boundary",fa$boundary,"LL",fa$log_likelihood,"\n")
d2 <- subset(d, g %in% which(sizes>1))
fb <- caught(multilpa(d2,c("a","b","c"),"g",n_profiles=3,n_group_classes=1,centering="person",n_starts=5,seed=1))
cat("without singletons: warnings:",attr(fb,"w"),"\n"); print(round(fb$means,3)); print(round(fb$variances,5)); print(round(drop(fb$profile_probabilities),3)); cat("boundary",fb$boundary,"LL",fb$log_likelihood,"\n")
fa2 <- caught(multilpa(d,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,centering="person",n_starts=5,seed=1))
fb2 <- caught(multilpa(d2,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,centering="person",n_starts=5,seed=1))
cat("k=2 means diff with/without singletons",max(abs(fa2$means[order(fa2$means[,1]),]-fb2$means[order(fb2$means[,1]),])),"warnings",attr(fa2,"w"),"\n")
print(round(fa2$variances,4)); print(round(fb2$variances,4)); print(round(fa2$means,4))
# centring with one-row-per-group only on one indicator through missingness: FIML
d3 <- d; d3$c[d3$g %in% 13:40 & duplicated(d3$g)] <- NA
fc <- caught(multilpa(d3,c("a","b","c"),"g",n_profiles=2,n_group_classes=1,centering="person",missing="fiml",n_starts=3,seed=1))
cat("indicator observed once in 40 groups: warnings:",attr(fc,"w")," zero values on c:",sum(fc$indicator_data[,"c"]==0,na.rm=TRUE),"\n")
# Same data, same number of profiles, spike excluded: how much likelihood the fabricated zeros buy
fx <- caught(multilpa(d,c("a","b","c"),"g",n_profiles=3,n_group_classes=1,centering="person",n_starts=5,seed=1,boundary="exclude"))
cat("k=3 same data: spike LL",fa$log_likelihood," excluding bound",fx$log_likelihood," gain",fa$log_likelihood-fx$log_likelihood,"\n")
cat("BIC k=2",fa2$bic," k=3",fa$bic," (true k=2); without singletons k=2",fb2$bic," k=3",fb$bic,"\n")
cat("rows with all-zero deviations:",sum(rowSums(abs(fa$indicator_data))==0),"; modal in spike profile:",sum(fa$subject_profiles[rowSums(abs(fa$indicator_data))==0]==which.min(fa$variances[,1])),"\n")
