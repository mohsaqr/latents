suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(21)
gen <- function() { nid <- 100; nt <- 3; st <- matrix(0L,nid,nt); st[,1] <- 1+rbinom(nid,1,.5); for (t in 2:nt) st[,t] <- ifelse(runif(nid)<.85,st[,t-1],3-st[,t-1])
  dl <- data.frame(id=rep(seq_len(nid),each=nt), time=rep(seq_len(nt),nid)); s <- as.vector(t(st)); dl$s <- s
  dl$y1 <- rnorm(nid*nt, c(-2,2)[s]); dl$y2 <- rnorm(nid*nt, c(0,3)[s]); dl }
dl <- gen()
hl <- lta(dl, c("y1","y2"), "id", n_profiles=2, time="time", model="VVI", n_starts=3, seed=1)
print(cor(dl$y1[dl$s==1],dl$y2[dl$s==1])); print(cor(dl$y1[dl$s==2],dl$y2[dl$s==2]))
print(hl$means); print(hl$variances); print(colSums(hl$subject_posteriors))
print(get_results(hl,"residuals")[,c("profile","observed","expected","effective_n","statistic","p_value")])
P <- hl$subject_posteriors; for (k in 1:2) { e1 <- dl$y1-hl$means[k,1]; e2 <- dl$y2-hl$means[k,2]; w <- P[,k]; cat("hand r", k, sum(w*e1*e2)/sqrt(sum(w*e1^2)*sum(w*e2^2)), "\n") }
# rerun several data sets, same generator; compare the plain LPA as a control
res <- t(sapply(1:25, function(i) { set.seed(100+i); d <- gen()
  f <- suppressWarnings(lta(d, c("y1","y2"), "id", n_profiles=2, time="time", model="VVI", n_starts=2, seed=1)); r <- get_results(f,"residuals")
  p <- suppressWarnings(lpa(d, c("y1","y2"), n_profiles=2, model="VVI", n_starts=2, seed=1)); q <- get_results(p,"residuals")
  c(lta = min(r$p_value), lpa = min(q$p_value), lta_abs = max(abs(r$statistic)), lpa_abs = max(abs(q$statistic))) }))
print(round(colMeans(res[,1:2] < .05),3)); print(round(colMeans(res[,3:4]),3))
