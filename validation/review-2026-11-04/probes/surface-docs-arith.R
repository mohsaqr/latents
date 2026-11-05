# Arithmetic behind the article's sentence on AvePP against entropy (base R only).
re <- function(P, w = rep(1, nrow(P))) { 1 - sum(w * -rowSums(P * log(pmax(P, 1e-300)))) / (sum(w) * log(ncol(P))) }
ave <- function(P, w = rep(1, nrow(P))) sum(w * apply(P, 1, max)) / sum(w)
tail2 <- rbind(c(.99,.01), c(.5,.5)); wt <- c(.9,.1)
cat(sprintf("K=2  90%% at .99, 10%% flat : relative entropy %.3f, mean max posterior %.3f\n", re(tail2, wt), ave(tail2, wt)))
tail4 <- rbind(c(.97,.01,.01,.01), c(.7,.1,.1,.1)); wt4 <- c(.8,.2)
cat(sprintf("K=4  80%% at .97, 20%% at .7 : relative entropy %.3f, mean max posterior %.3f\n", re(tail4, wt4), ave(tail4, wt4)))
flat <- rbind(c(.9,.1)); cat(sprintf("K=2  all at .90           : relative entropy %.3f, mean max posterior %.3f\n", re(flat), ave(flat)))
flat <- rbind(c(.8,.2)); cat(sprintf("K=2  all at .80           : relative entropy %.3f, mean max posterior %.3f\n", re(flat), ave(flat)))
cat(sprintf("SABIC per-parameter penalty at n = 106: %.3f; at n = 1422: %.3f\n", log(108/24), log(1424/24)))
cat(sprintf("Staged BIC with 23 parameters: %.2f\n", 2*8440.08 + 23*log(106)))
