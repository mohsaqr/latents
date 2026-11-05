suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(3)
ng <- 50; nper <- 6
g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each = nper)
k <- 1 + rbinom(length(id), 1, ifelse(g[id] == 1, .8, .2))
d <- data.frame(id = id, y1 = rnorm(length(id), c(-1, 1)[k]), y2 = rnorm(length(id), c(0, 1.5)[k]), y3 = rnorm(length(id)))
d$a <- ifelse(runif(nrow(d)) < c(.2,.8)[k], "x", "y"); d$b <- ifelse(runif(nrow(d)) < c(.8,.2)[k], "u", "v"); d$c <- ifelse(runif(nrow(d)) < c(.3,.7)[k], "r", "s")
f <- multilpa(d, c("y1","y2","y3"), "id", n_profiles = 2, n_group_classes = 2, model = "VVI", n_starts = 5, seed = 1)
# --- classification table by hand
P <- f$subject_posteriors; G <- f$group_posteriors
cl <- get_results(f, "classification")
hand <- function(P) { a <- max.col(P, "first"); K <- ncol(P)
  avepp <- vapply(1:K, function(j) mean(P[a == j, j]), 0); pi <- colMeans(P)
  data.frame(n = tabulate(a, K), est_n = colSums(P), avepp = avepp, occ = (avepp/(1-avepp))/(pi/(1-pi))) }
print(cl); print(hand(P)); print(hand(G))
cat("subject_profiles equals max.col:", identical(as.integer(f$subject_profiles), max.col(P, "first")), " group:", identical(as.integer(f$group_classes), max.col(G, "first")), "\n")
# --- error matrix and BCH by hand
errs <- get_results(f, "classification_errors"); print(errs)
a <- max.col(P, "first"); D <- matrix(0, 2, 2)
for (t in 1:2) for (s in 1:2) D[t, s] <- sum(P[, t] * (a == s)) / sum(P[, t])
print(D)
b <- get_results(f, "bch_weights", level = "individuals")
Dinv <- solve(D)
hb <- as.vector(t(Dinv[a, ]))
cat("BCH max diff to hand:", max(abs(b$weight - hb)), "\n")
# BCH logic: on the modal table, weights should reproduce posterior class total
#   sum_i w_it = sum_s n_s Dinv[s,t]; property: sum_i w_it = n * pi_t when n_s = sum_t N_t D_ts
Nt <- colSums(P); ns <- tabulate(a, 2)
cat("n_s - t(D)%*%N_t:", as.vector(ns - t(D) %*% Nt), "\n")
cat("sum BCH weights by class:", tapply(b$weight, b$class, sum), " posterior totals:", Nt, "\n")
# Independent BCH unbiasedness by simulation? small: one dataset with known outcome effect
# --- bivariate residuals by hand
r <- get_results(f, "residuals"); print(head(r, 4))
ro <- r[r$profile == "profile_1" & r$indicator_1 == "y1" & r$indicator_2 == "y2", ]
w <- P[, 1]; e1 <- d$y1 - f$means[1,1]; e2 <- d$y2 - f$means[1,2]
ob <- sum(w*e1*e2)/sqrt(sum(w*e1^2)*sum(w*e2^2)); z <- atanh(ob)*sqrt(sum(w)-3)
cat("gaussian pair: obs", ro$observed, "hand", ob, " stat", ro$statistic, "hand", z, " p", ro$p_value, "hand", 2*pnorm(-abs(z)), "\n")
lc <- multilpa(d, c("a","b","c"), "id", n_profiles = 2, n_group_classes = 2, categorical = c("a","b","c"), n_starts = 5, seed = 1)
rc <- get_results(lc, "residuals"); print(head(rc, 3))
Pc <- lc$subject_posteriors
tab <- function(w, A, B) { o <- tapply(w, list(A, B), sum); o[is.na(o)] <- 0; o }
o <- tab(Pc[,1], d$a, d$b); pa <- lc$response_probabilities$a[1,]; pb <- lc$response_probabilities$b[1,]
ex <- sum(Pc[,1]) * outer(pa[rownames(o)], pb[colnames(o)])
x2 <- sum((o-ex)^2/ex); rc1 <- rc[rc$profile=="profile_1" & rc$indicator_1=="a" & rc$indicator_2=="b",]
cat("cat pair: stat", rc1$statistic, "hand", x2, " df", rc1$df, " p", rc1$p_value, "hand", pchisq(x2,1,lower.tail=FALSE), "\n")
rc2 <- get_results(lc, "residuals", by = "overall"); print(rc2[rc2$indicator_1=="a" & rc2$indicator_2=="b",])
PP <- colSums(Pc); ov <- table(d$a, d$b)
ex2 <- PP[1]*outer(lc$response_probabilities$a[1,], lc$response_probabilities$b[1,]) + PP[2]*outer(lc$response_probabilities$a[2,], lc$response_probabilities$b[2,])
cat("overall hand", sum((ov-ex2[rownames(ov),colnames(ov)])^2/ex2[rownames(ov),colnames(ov)]), "\n")
# wording: "largest" vs lowest p
qd <- diagnostics(f, plots = FALSE); ck <- get_results(qd, "checks")
print(ck[ck$domain == "Local independence", c("check","display","criterion","status","detail")])
rd <- get_results(qd, "residuals"); print(rd[1:3,c("profile","indicator_1","indicator_2","residual","effective_n","statistic","p_value","p_adjusted")]); print(rd[which.min(rd$p_adjusted),c("profile","indicator_1","indicator_2","residual","effective_n","p_adjusted")])
print(ck[,c("domain","check","level","display","status")])
# BIC row
cat("info BIC", f$bic, "is groups:", isTRUE(all.equal(f$bic, get_results(f,"information_criteria")$bic_groups)), "\n")
print(compare_models(a = f, b = multilpa(d, c("y1","y2","y3"), "id", n_profiles = 2, n_group_classes = 1, model = "VVI", n_starts = 5, seed = 1)))
# --- tree flows
gridw <- enumerate_classes(d, c("y1","y2"), "id", n_profiles = 1:3, n_group_classes = 1:2, model = "VVI", n_starts = 3, seed = 1)
tr <- get_results(gridw, "tree"); print(head(tr))
chk <- do.call(rbind, by(tr, list(tr$tree, tr$from_candidate, tr$to_candidate), function(s) {
  fm <- tapply(s$mass, s$from_class, sum); tm <- tapply(s$mass, s$to_class, sum)
  data.frame(from = s$from_candidate[1], to = s$to_candidate[1],
    total = sum(s$mass), from_gap = max(abs(fm - s$from_share[match(names(fm), s$from_class)])), to_gap = max(abs(tm - s$to_share[match(names(tm), s$to_class)])))
}))
print(chk)
# hand share against fit
fit3 <- get_candidate(gridw, n_profiles = 3, n_group_classes = 2); fit2 <- get_candidate(gridw, n_profiles = 2, n_group_classes = 2)
M <- crossprod(fit2$subject_posteriors, fit3$subject_posteriors)/nrow(d); print(round(M,4))
print(tr[tr$from_candidate == .latents_candidate_names(gridw)[gridw$table$n_profiles==2 & gridw$table$n_group_classes==2] & tr$to_candidate == .latents_candidate_names(gridw)[gridw$table$n_profiles==3 & gridw$table$n_group_classes==2], c("from_class","to_class","mass","share_of_from","share_of_to")])
gr_tree <- get_results(gridw, "tree", level = "groups"); 
chk <- do.call(rbind, by(gr_tree, list(gr_tree$tree, gr_tree$from_candidate, gr_tree$to_candidate), function(s) data.frame(total=sum(s$mass), from = s$from_candidate[1], to = s$to_candidate[1])))
print(chk)
# weighted and FIML tree
d$w <- rep(sample(c(0,1,2,3), ng, TRUE, prob = c(.1,.3,.3,.3)), each = nper)
dm <- d; dm$y1[sample(nrow(d), 30)] <- NA; dm$y2[which(is.na(dm$y1))[1:8]] <- NA
gw <- enumerate_classes(dm, c("y1","y2"), "id", n_profiles = 1:3, n_group_classes = 2, model = "VVI", n_starts = 3, seed = 1, weights = "w", missing = "fiml")
tw <- get_results(gw, "tree"); 
chk <- do.call(rbind, by(tw, list(tw$from_candidate, tw$to_candidate), function(s) {
  fm <- tapply(s$mass, s$from_class, sum); data.frame(from = s$from_candidate[1], to = s$to_candidate[1], total = sum(s$mass), gap = max(abs(fm - s$from_share[match(names(fm), s$from_class)])))
})); print(chk)
fw <- get_candidate(gw, n_profiles = 2, n_group_classes = 2)
print(get_results(fw, "classification")); print(get_results(fw, "entropy"))
cat("n_informative", fw$n_informative, "n_informative_groups", fw$n_informative_groups, " rows with w>0:", sum(dm$w > 0), " rows w>0 and observed:", sum(dm$w > 0 & !(is.na(dm$y1) & is.na(dm$y2))), " groups w>0:", length(unique(dm$id[dm$w>0])), "\n")
wi <- dm$w; wn <- fw$sampling_weights[fw$group_index]; keep <- wi > 0 & !(is.na(dm$y1) & is.na(dm$y2))
Pw <- fw$subject_posteriors[keep,]; cat("hand relent:", 1 - sum(wn[keep] * -rowSums(ifelse(Pw>0,Pw*log(Pw),0)))/(sum(wn[keep])*log(2)), "\n")
print(get_results(fw, "information_criteria", format="long")[c(4,10),])
# modal-vs-estimated gap under weights (separated classes; weight tied to class)
set.seed(4); n <- 400; kk <- rep(1:2, c(n/2,n/2)); ds <- data.frame(id = seq_len(n), y1 = rnorm(n, c(-6,6)[kk]), y2 = rnorm(n, c(-6,6)[kk]), w = c(1,5)[kk])
sw <- lpa(ds, c("y1","y2"), n_profiles = 2, model = "EEI", weights = "w", n_starts = 2, seed = 1)
print(get_results(sw, "classification")); ckw <- get_results(diagnostics(sw, plots = FALSE), "checks"); print(ckw[ckw$check %in% c("modal vs estimated proportions","smallest class","relative entropy"), c("check","display","status")])
# boundary flag on probabilities
set.seed(5); n <- 300; kk <- rbinom(n,1,.5)+1
dp <- data.frame(a = ifelse(kk==1, "x", ifelse(runif(n)<.5,"x","y")), b = ifelse(kk==1,"u", ifelse(runif(n)<.3,"u","v")), c = ifelse(kk==1, ifelse(runif(n)<.8,"r","s"), ifelse(runif(n)<.2,"r","s")), d = ifelse(kk==1, ifelse(runif(n)<.9,"p","q"), ifelse(runif(n)<.3,"p","q")))
lp <- lca(dp, c("a","b","c","d"), n_classes = 2, n_starts = 5, seed = 1)
print(lp$response_probabilities$a); print(lp$response_probabilities$b)
cat("boundary:", lp$boundary, "\n")
ckl <- get_results(diagnostics(lp, plots = FALSE), "checks"); print(ckl[ckl$check == "estimates at a bound", c("check","display","status","detail")])
pi0 <- tryCatch(parameter_inference(lp), error = function(e) e); print(class(pi0)[1]); if (inherits(pi0,"error")) cat(conditionMessage(pi0),"\n")
print(head(tryCatch(get_results(lp,"responses"), error=function(e) conditionMessage(e)),4))
