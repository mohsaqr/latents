suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(8)
ng <- 40; nper <- 5; g <- sample(1:2, ng, TRUE); id <- rep(seq_len(ng), each=nper); k <- 1+rbinom(length(id),1,ifelse(g[id]==1,.8,.2))
d <- data.frame(id=id, y1=rnorm(length(id),c(-1.5,1.5)[k]))
d$a <- ifelse(runif(nrow(d)) < c(.2,.8)[k], "x", "y"); d$b <- ifelse(runif(nrow(d)) < c(.7,.3)[k], "u", "v")
s <- list(vars=c("a","b","y1"), categorical=c("a","b"), local_dependence="a:b", model="VVI")
mk <- function(K, G) do.call(multilpa, c(list(data=d, id="id", n_profiles=K, n_group_classes=G, n_starts=2, seed=1), s))
n <- mk(2,1); alt <- mk(2,2)
str(n[c("vars","original_vars","categorical")])
print(names(n$response_probabilities))
r <- tryCatch(bootstrap_lrt(n, alt, iter=3, seed=1), error=function(e) e)
print(class(r)); if (inherits(r,"error")) cat(conditionMessage(r), "\n") else { print(as.data.frame(r)); print(get_results(r,"replicates")) }
cat("---- simulate direct\n")
sim <- tryCatch(.multilpa_simulate(n), error=function(e) e); print(class(sim)); if (!inherits(sim,"error")) print(head(sim,3))
# simulated margins and the moment checks for NB
# non-syntactic categorical name, no joint item
d2 <- d; names(d2)[names(d2) == "a"] <- "item one"; names(d2)[names(d2) == "b"] <- "item-two"
mk2 <- function(K) lca(d2, c("item one","item-two"), n_classes = K, n_starts = 2, seed = 1)
s1 <- mk2(1); s2 <- mk2(2)
r2 <- tryCatch(bootstrap_lrt(s1, s2, iter = 3, seed = 1), error = function(e) e, warning = function(w) w)
print(class(r2)); cat(conditionMessage(r2), "\n")
r2 <- suppressWarnings(tryCatch(bootstrap_lrt(s1, s2, iter = 3, seed = 1), error = function(e) e)); if (!inherits(r2,"error")) { print(as.data.frame(r2)[,c("statistic","p_value","n_valid","n_invalid")]); print(get_results(r2,"replicates")[1,c("statistic","valid","error")]) }
# continuous non-syntactic names
d3 <- d; names(d3)[names(d3) == "y1"] <- "score 1"; d3$`score 2` <- rnorm(nrow(d3))
c1 <- lpa(d3, c("score 1","score 2"), n_profiles = 1, model = "EEI", n_starts = 2, seed = 1); c2 <- lpa(d3, c("score 1","score 2"), n_profiles = 2, model = "EEI", n_starts = 2, seed = 1)
r3 <- suppressWarnings(tryCatch(bootstrap_lrt(c1, c2, iter = 3, seed = 1), error = function(e) e)); if (!inherits(r3,"error")) print(as.data.frame(r3)[,c("statistic","p_value","n_valid","n_invalid")]) else cat(conditionMessage(r3),"\n")
d$cnt <- rnbinom(nrow(d), size=2, mu=c(1,5)[k])
nb <- multilpa(d, c("cnt","y1"), "id", n_profiles=2, n_group_classes=2, count="cnt", count_model="negative_binomial", count_dispersion="equal", model="VVI", n_starts=2, seed=1)
sim <- do.call(rbind, lapply(1:50, function(i) .multilpa_simulate(nb)))
wk <- as.vector(nb$group_probabilities %*% nb$profile_probabilities); mu <- nb$count_means[,1]; disp <- nb$count_dispersion[1,1]
cat("NB sim mean", mean(sim$cnt), "model", sum(wk*mu), " var sim", var(sim$cnt), "model", sum(wk*(mu+disp*mu^2+mu^2))-sum(wk*mu)^2, "\n")
# ---- deltas, knees, Schwarz by hand
dv2 <- data.frame(id = rep(1:30, each = 4)); kk <- rep(rbinom(30,1,.5), each=4); dv2$y1 <- rnorm(120, 2*kk); dv2$y2 <- rnorm(120, -2*kk)
gg <- enumerate_lpa(dv2, c("y1","y2"), n_profiles = 1:5, model = "EEI", n_starts = 3, seed = 1)
tb <- get_results(gg, "candidates"); print(tb[,c("candidate","bic_groups","bic_individual")])
dl <- get_results(gg, "deltas"); print(dl[,c("candidate","value","delta_profile","from_profile")])
kn <- get_results(gg, "knees"); print(kn[,c("direction","knee_candidate","distance","knee_without_largest","stable")])
x <- 1:5; y <- tb$bic_groups; xs <- (x-1)/4; ys <- (y-min(y))/(max(y)-min(y)); sl <- ys[5]-ys[1]; dist <- ((ys[1]+sl*xs)-ys)/sqrt(1+sl^2); cat("hand distances", round(dist,4), "\n")
cm <- as.data.frame(compare_models(gg)); dlt <- tb$bic_groups-min(tb$bic_groups); wgt <- exp(-dlt/2)/sum(exp(-dlt/2)); cat("Schwarz diff", max(abs(cm$bic_weight-wgt)), "\n")
print(get_results(gg,"best")[,c("criterion","convention","candidate","value")])
