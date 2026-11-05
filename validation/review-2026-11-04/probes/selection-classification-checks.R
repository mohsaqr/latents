suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(11)
# (a) LCA with a probability on its bound: fit$boundary, diagnostics, grid, LRT, BLRT, LMR
n <- 400; kk <- rbinom(n,1,.5)+1
dp <- data.frame(id = seq_len(n), a = ifelse(kk==1, "x", ifelse(runif(n)<.4,"x","y")), b = ifelse(runif(n)<c(.9,.2)[kk],"u","v"), c = ifelse(runif(n)<c(.85,.2)[kk],"r","s"), d = ifelse(runif(n)<c(.9,.3)[kk],"p","q"), e = ifelse(runif(n)<c(.15,.7)[kk],"m","n"))
v <- c("a","b","c","d","e")
lp <- lca(dp, v, n_classes = 2, n_starts = 8, seed = 1)
cat("response a by class:\n"); print(lp$response_probabilities$a)
cat("min_probability", lp$min_probability, " fit$boundary", lp$boundary, "\n")
ck <- get_results(diagnostics(lp, plots = FALSE), "checks"); print(ck[ck$check == "estimates at a bound", c("check","display","status")])
pi0 <- tryCatch(parameter_inference(lp), error = function(e) e); cat("parameter_inference:", class(pi0)[1], substr(if (inherits(pi0,"error")) conditionMessage(pi0) else "", 1, 100), "\n")
lj <- lca(dp, v, n_classes = 2, n_starts = 8, seed = 1, local_dependence = "b:c")
r <- withCallingHandlers(as.data.frame(compare_models(free = lp, joint = lj, test = "lrt"))[,c("model","n_parameters","log_likelihood","lr_statistic","lr_df","lr_p_value")],
  warning = function(w) { cat("WARNING:", class(w)[1], conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(r); cat("joint fit boundary:", lj$boundary, " any prob at bound:", any(unlist(lj$response_probabilities) <= lj$min_probability*(1+1e-8)), "\n")
gr <- enumerate_lca(dp, v, n_classes = 1:3, n_starts = 8, seed = 1)
print(get_results(gr,"candidates")[,c("candidate","boundary","converged")])
lm <- get_results(gr, "lmr"); print(lm[,c("null_candidate","alternative_candidate","vlmr_p","almr_p","note")])
# NB boundary flag
set.seed(2); dn <- data.frame(id = rep(1:40, each=5)); kn <- 1+rbinom(200,1,.5); dn$cnt <- rpois(200, c(1,5)[kn]); dn$y <- rnorm(200, c(-1,1)[kn])
nb <- suppressWarnings(multilpa(dn, c("cnt","y"), "id", n_profiles=2, n_group_classes=1, count="cnt", count_model="negative_binomial", model="VVI", n_starts=3, seed=1))
cat("NB dispersion:", nb$count_dispersion, " boundary flag", nb$boundary, "\n"); ckn <- get_results(diagnostics(nb, plots=FALSE),"checks"); print(ckn[ckn$check=="estimates at a bound",c("display","status")])
# Homogeneous LTA probability bound
set.seed(3); nid <- 80; nt <- 3; st <- matrix(0L,nid,nt); st[,1] <- 1+rbinom(nid,1,.5); for (t in 2:nt) st[,t] <- ifelse(runif(nid)<.85,st[,t-1],3-st[,t-1])
dl <- data.frame(id=rep(seq_len(nid),each=nt), time=rep(seq_len(nt),nid)); s <- as.vector(t(st))
dl$a <- ifelse(s==1,"x",ifelse(runif(nid*nt)<.4,"x","y")); dl$b <- ifelse(runif(nid*nt)<c(.9,.2)[s],"u","v"); dl$c <- ifelse(runif(nid*nt)<c(.85,.2)[s],"r","s")
hl <- lta(dl, c("a","b","c"), "id", n_profiles=2, time="time", categorical=c("a","b","c"), n_starts=5, seed=1)
cat("LTA class", class(hl), "; a probs", round(as.vector(hl$response_probabilities$a),6), " boundary", hl$boundary, "\n")
ol <- lta(dl, c("a","b","c"), "id", n_profiles=2, time="time", categorical=c("a","b","c"), transitions="occasion", n_starts=5, seed=1)
cat("General LTA class", class(ol), " boundary", ol$boundary, "\n")
# (b) drop_empty FALSE sample size
set.seed(4); dm <- data.frame(id=rep(1:30,each=4), y1=rnorm(120), y2=rnorm(120)); dm$y1[1:4] <- NA; dm$y2[1:4] <- NA; dm$y1[9] <- NA; dm$y2[9] <- NA
fD <- multilpa(dm, c("y1","y2"), "id", n_profiles=2, n_group_classes=1, model="EEI", missing="fiml", n_starts=2, seed=1)
fK <- multilpa(dm, c("y1","y2"), "id", n_profiles=2, n_group_classes=1, model="EEI", missing="fiml", n_starts=2, seed=1, drop_empty=FALSE)
print(rbind(drop=get_results(fD,"information_criteria",format="long")[c(4,10),"n"], keep=get_results(fK,"information_criteria",format="long")[c(4,10),"n"]))
cat("n_informative", fD$n_informative, fK$n_informative, "; rows 120, groups 30\n")
print(rbind(get_results(fD,"entropy"), get_results(fK,"entropy")))
# (c) largest wording: small profile with large residual
set.seed(5); nn <- 600; kk <- ifelse(runif(nn)<.1,1,2); ds <- data.frame(id=seq_len(nn), y1=rnorm(nn,c(-5,3)[kk]), y2=rnorm(nn,c(-5,3)[kk]), y3=rnorm(nn))
e <- rnorm(nn); ds$y2[kk==2] <- ds$y2[kk==2] + .35*(ds$y1[kk==2]-3)
sf <- lpa(ds, c("y1","y2","y3"), n_profiles=2, model="VVI", n_starts=3, seed=1)
qd <- diagnostics(sf, plots=FALSE); ck <- get_results(qd,"checks"); cat(ck$detail[ck$domain=="Local independence"], "\n")
print(get_results(qd,"residuals")[1:4,c("profile","indicator_1","indicator_2","residual","effective_n","p_adjusted")])
# (d) general LTA with group classes in diagnostics
set.seed(6); ng <- 40; dg <- do.call(rbind, lapply(seq_len(ng), function(j) { h <- sample(1:2,1); do.call(rbind, lapply(1:6, function(p) { st <- 1+rbinom(3,1,c(.2,.8)[h]); data.frame(cl=j, id=j*100+p, time=1:3, y1=rnorm(3,c(-2,2)[st]), y2=rnorm(3,c(0,3)[st])) })) }))
lg <- lta(dg, c("y1","y2"), "id", n_profiles=2, time="time", n_group_classes=2, transitions="occasion", model="VVI", n_starts=3, seed=1)
cat("general LTA n_group_classes field:", format(lg$n_group_classes), " arg", lg$arguments$n_group_classes, " posteriors", ncol(lg$group_posteriors), "\n")
ck <- get_results(diagnostics(lg, plots=FALSE),"checks"); print(ck[,c("check","level","display","status")])
print(tail(get_results(lg,"fit"),1))
ldf <- compare_models(a=lg, b=lta(dg, c("y1","y2"), "id", n_profiles=2, time="time", n_group_classes=2, transitions="occasion", model="EEI", n_starts=3, seed=1)); print(as.data.frame(ldf)[,c("model","n_group_classes","n_parameters","bic_groups","icl_groups","icl_individual","profile_entropy","group_entropy")])
# (e) joint-pair BLRT with original data
set.seed(8); ng <- 40; nper <- 5; g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each=nper); k <- 1+rbinom(length(id),1,ifelse(g[id]==1,.8,.2))
d <- data.frame(id=id, y1=rnorm(length(id),c(-1.5,1.5)[k])); d$a <- ifelse(runif(nrow(d)) < c(.2,.8)[k], "x", "y"); d$b <- ifelse(runif(nrow(d)) < c(.7,.3)[k], "u", "v")
s <- list(vars=c("a","b","y1"), categorical=c("a","b"), local_dependence="a:b", model="VVI")
mk <- function(K, G) do.call(multilpa, c(list(data=d, id="id", n_profiles=K, n_group_classes=G, n_starts=2, seed=1), s))
n0 <- mk(2,1); a0 <- mk(2,2)
r <- tryCatch(bootstrap_lrt(n0, a0, data = d, iter=2, seed=1), error=function(e) e); cat("BLRT with data=:", class(r)[1], if (inherits(r,"error")) conditionMessage(r), "\n")
# report text
f <- multilpa(d, c("y1"), "id", n_profiles=2, n_group_classes=2, n_starts=3, seed=1, model = "VVI")
report(f, plots=FALSE)
