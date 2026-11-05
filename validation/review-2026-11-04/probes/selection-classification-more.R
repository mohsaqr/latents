suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(21)
ns <- asNamespace("latents")
# (1) general and homogeneous LTA residuals: table and the check's reason
nid <- 100; nt <- 3; st <- matrix(0L,nid,nt); st[,1] <- 1+rbinom(nid,1,.5); for (t in 2:nt) st[,t] <- ifelse(runif(nid)<.85,st[,t-1],3-st[,t-1])
dl <- data.frame(id=rep(seq_len(nid),each=nt), time=rep(seq_len(nt),nid)); s <- as.vector(t(st))
dl$y1 <- rnorm(nid*nt, c(-2,2)[s]); dl$y2 <- rnorm(nid*nt, c(0,3)[s])
hl <- lta(dl, c("y1","y2"), "id", n_profiles=2, time="time", model="VVI", n_starts=3, seed=1)
ol <- lta(dl, c("y1","y2"), "id", n_profiles=2, time="time", model="VVI", transitions="occasion", n_starts=3, seed=1)
for (f in list(hl, ol)) {
  cat("----", class(f), "\n")
  r <- tryCatch(.multilpa_bivariate_residuals(f, dl), error = function(e) e)
  if (inherits(r,"error")) cat("internal residuals:", class(r)[1], conditionMessage(r), "\n") else cat("internal residual rows:", nrow(r), "\n")
  print(class(tryCatch(get_results(f, "residuals"), error = function(e) e))[1])
  q <- diagnostics(f, data = dl, plots = FALSE); ck <- get_results(q,"checks"); print(ck[ck$domain=="Local independence",c("check","display","status","detail")])
  rr <- tryCatch(get_residuals(f), error=function(e) e); cat("get_residuals:", class(rr)[1], if (inherits(rr,"error")) conditionMessage(rr) else nrow(rr), "\n")
}
# (2) single-indicator profile fit: the reason text
d1 <- data.frame(id=rep(1:40,each=5), y=rnorm(200))
f1 <- multilpa(d1, "y", "id", n_profiles=2, n_group_classes=2, n_starts=2, seed=1)
ck <- get_results(diagnostics(f1, plots=FALSE),"checks"); print(ck[ck$domain=="Local independence",c("check","status","detail")])
# (3) wording: "largest" names the smallest p, not the largest residual
set.seed(5); nn <- 800; kk <- ifelse(runif(nn)<.08,1,2); ds <- data.frame(id=seq_len(nn), y1=rnorm(nn,c(-6,4)[kk]), y2=rnorm(nn,c(-6,4)[kk]))
big <- which(kk == 2); small <- which(kk == 1)
ds$y2[big] <- ds$y2[big] + .22*(ds$y1[big]-4); ds$y2[small] <- ds$y2[small] + .45*(ds$y1[small]+6)
sf <- lpa(ds, c("y1","y2"), n_profiles=2, model="VVI", n_starts=3, seed=1)
qd <- diagnostics(sf, plots=FALSE); ck <- get_results(qd,"checks"); cat(ck$detail[ck$domain=="Local independence"], "\n")
print(get_results(qd,"residuals")[,c("profile","indicator_1","indicator_2","residual","effective_n","p_adjusted")])
# (4) fix for draw: check.names = FALSE, in memory only
draw_src <- body(ns$.multilpa_draw_indicators)
txt <- deparse(ns$.multilpa_draw_indicators)
txt2 <- sub("as.data.frame(drawn, stringsAsFactors = FALSE)", "as.data.frame(drawn, stringsAsFactors = FALSE, check.names = FALSE)", txt, fixed = TRUE)
stopifnot(!identical(txt, txt2))
fn <- eval(parse(text = txt2), envir = ns); environment(fn) <- ns
unlockBinding(".multilpa_draw_indicators", ns); assign(".multilpa_draw_indicators", fn, envir = ns); lockBinding(".multilpa_draw_indicators", ns)
set.seed(8); ng <- 40; nper <- 5; g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each=nper); k <- 1+rbinom(length(id),1,ifelse(g[id]==1,.8,.2))
d <- data.frame(id=id, y1=rnorm(length(id),c(-1.5,1.5)[k])); d$a <- ifelse(runif(nrow(d)) < c(.2,.8)[k], "x", "y"); d$b <- ifelse(runif(nrow(d)) < c(.7,.3)[k], "u", "v")
s <- list(vars=c("a","b","y1"), categorical=c("a","b"), local_dependence="a:b", model="VVI")
mk <- function(K, G, data = d, ...) do.call(multilpa, c(list(data=data, id="id", n_profiles=K, n_group_classes=G, n_starts=2, seed=1), s, list(...)))
n0 <- mk(2,1); a0 <- mk(2,2)
r <- bootstrap_lrt(n0, a0, iter=5, seed=1); print(as.data.frame(r)[,c("statistic","p_value","n_valid","n_invalid")]); print(get_results(r,"replicates")[,c("statistic","valid","error")])
sim <- do.call(rbind, lapply(1:40, function(i) .multilpa_simulate(n0)))
print(round(prop.table(table(sim[["a:b"]])),3)); wk <- as.vector(n0$group_probabilities %*% n0$profile_probabilities); print(round(colSums(wk * n0$response_probabilities[["a:b"]]),3))
# FIML joint, half observed pair after the fix
dm <- d; dm$a[1:6] <- NA; dm$b[7:9] <- NA; dm$y1[10:14] <- NA
nm <- mk(2,1,data=dm,missing="fiml"); am <- mk(2,2,data=dm,missing="fiml")
rm <- tryCatch(bootstrap_lrt(nm, am, iter=3, seed=1), error=function(e) e, warning=function(w) w)
print(class(rm)[1]); if (inherits(rm,"condition")) cat(conditionMessage(rm),"\n") else print(as.data.frame(rm)[,c("statistic","p_value","n_valid")])
rm <- suppressWarnings(tryCatch(bootstrap_lrt(nm, am, iter=3, seed=1), error=function(e) e)); if (!inherits(rm,"error")) print(get_results(rm,"replicates")[,c("statistic","valid","error")])
# (5) One data set: LCA BVR power with real within-class dependence
set.seed(31); n <- 800; kk <- rbinom(n,1,.5)+1
# a and b dependent within class 1 (phi about .4), independent otherwise
la <- rbinom(n,1,c(.5,.8)[kk]); lb <- ifelse(kk==1, ifelse(runif(n)<.7, la, rbinom(n,1,.5)), rbinom(n,1,.2))
dc <- data.frame(a=la, b=lb, c=rbinom(n,1,c(.8,.2)[kk]), d=rbinom(n,1,c(.75,.25)[kk]), e=rbinom(n,1,c(.2,.8)[kk]))
fc <- lca(dc, c("a","b","c","d","e"), n_classes=2, n_starts=5, seed=1)
print(cor(dc$a[kk==1], dc$b[kk==1]))
print(get_results(fc,"residuals")[1:4,c("profile","indicator_1","indicator_2","observed","statistic","df","p_value")])
print(get_results(diagnostics(fc, plots=FALSE),"checks")[,c("check","display","status")][16,])
ckk <- get_results(diagnostics(fc, plots=FALSE),"checks"); print(ckk[ckk$domain=="Local independence",c("display","status","detail")])
print(get_results(diagnostics(fc, plots=FALSE, by="overall"),"checks")[ckk$domain=="Local independence",c("display","status","detail")])
# (6) LRT silently chi-squares LCA with a probability on its bound
set.seed(9); n <- 500; kk <- rbinom(n,1,.5)+1
dq <- data.frame(id=seq_len(n), a=ifelse(kk==1,"x",ifelse(runif(n)<.5,"x","y")), b=ifelse(kk==1,"u",ifelse(runif(n)<.3,"u","v")), c=ifelse(runif(n)<c(.8,.2)[kk],"r","s"), d=ifelse(runif(n)<c(.9,.3)[kk],"p","q"))
l0 <- lca(dq, c("a","b","c","d"), n_classes=2, n_starts=6, seed=1)
l1 <- lca(dq, c("a","b","c","d"), n_classes=2, n_starts=6, seed=1, local_dependence="a:b")
cat("null at bound:", any(unlist(l0$response_probabilities) <= 1.0000001*l0$min_probability), " alt at bound:", any(unlist(l1$response_probabilities) <= 1.0000001*l1$min_probability), " flags", l0$boundary, l1$boundary, "\n")
print(l1$response_probabilities[["a:b"]])
tb <- withCallingHandlers(as.data.frame(compare_models(free=l0, joint=l1, test="lrt"))[,c("model","n_parameters","lr_statistic","lr_df","lr_p_value")], warning=function(w){cat("WARN:",class(w)[1],"\n");invokeRestart("muffleWarning")}); print(tb)
