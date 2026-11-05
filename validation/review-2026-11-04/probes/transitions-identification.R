# Does n_parameters count what the latent path distribution can identify?
# Pure mathematics: rank of d log p(path) / d (free logits) at random parameter
# values, compared with the count of free logits lta() uses. No package code.
set.seed(1)
softmax <- function(e) { e <- e - max(e); exp(e) / sum(exp(e)) }
path_probs <- function(theta, K, H, T, order, stayer, occasion) {
  pos <- 0L; take <- function(n) { v <- theta[pos + seq_len(n)]; pos <<- pos + n; v }
  shares <- softmax(c(0, take(H + stayer - 1L)))
  init <- lapply(seq_len(H + stayer), function(h) softmax(c(take(K - 1L), 0)))
  n_first <- if (occasion && order == 1L) T - 1L else 1L
  n_second <- if (order == 2L) (if (occasion) T - 2L else 1L) else 0L
  tr <- lapply(seq_len(H), function(h) lapply(seq_len(n_first), function(m)
    lapply(seq_len(K), function(k) { e <- numeric(K); e[-k] <- take(K - 1L); softmax(e) })))
  t2 <- lapply(seq_len(H), function(h) lapply(seq_len(n_second), function(m)
    lapply(seq_len(K^2), function(p) { j <- (p - 1L) %% K + 1L; e <- numeric(K); e[-j] <- take(K - 1L); softmax(e) })))
  paths <- as.matrix(expand.grid(rep(list(seq_len(K)), T)))
  apply(paths, 1L, function(path) {
    sum(vapply(seq_len(H + stayer), function(h) {
      p <- shares[h] * init[[h]][path[1]]
      for (t in 2:T) {
        if (h > H) { p <- p * (path[t] == path[t - 1L]); next }
        if (order == 2L && t >= 3L) {
          pair <- (path[t - 2L] - 1L) * K + path[t - 1L]
          p <- p * t2[[h]][[if (occasion) t - 2L else 1L]][[pair]][path[t]]
        } else p <- p * tr[[h]][[if (occasion && order == 1L) t - 1L else 1L]][[path[t - 1L]]][path[t]]
      }
      p
    }, 0))
  })
}
count <- function(K, H, T, order, stayer, occasion) {
  n_first <- if (occasion && order == 1L) T - 1L else 1L
  n_second <- if (order == 2L) (if (occasion) T - 2L else 1L) else 0L
  (H + stayer - 1L) + (H + stayer) * (K - 1L) + H * (n_first * K * (K - 1L) + n_second * K^2 * (K - 1L))
}
rank_at <- function(K, H, T, order, stayer, occasion) {
  n <- count(K, H, T, order, stayer, occasion)
  ranks <- vapply(1:3, function(s) {
    theta <- rnorm(n, 0, .7)
    jac <- vapply(seq_len(n), function(j) {
      e <- replace(numeric(n), j, 1e-5)
      (path_probs(theta + e, K, H, T, order, stayer, occasion) - path_probs(theta - e, K, H, T, order, stayer, occasion)) / 2e-5
    }, numeric(K^T))
    sv <- svd(jac)$d
    sum(sv > 1e-6 * sv[1])
  }, 0)
  max(ranks)
}
grid <- rbind(
  data.frame(K=2,H=1,T=2,order=1,stayer=0,occasion=FALSE),
  data.frame(K=2,H=2,T=2,order=1,stayer=0,occasion=FALSE),
  data.frame(K=3,H=2,T=2,order=1,stayer=0,occasion=FALSE),
  data.frame(K=2,H=2,T=3,order=1,stayer=0,occasion=FALSE),
  data.frame(K=3,H=2,T=3,order=1,stayer=0,occasion=FALSE),
  data.frame(K=2,H=2,T=4,order=1,stayer=0,occasion=FALSE),
  data.frame(K=2,H=3,T=4,order=1,stayer=0,occasion=FALSE),
  data.frame(K=2,H=2,T=5,order=1,stayer=0,occasion=FALSE),
  data.frame(K=2,H=1,T=2,order=1,stayer=1,occasion=FALSE),
  data.frame(K=3,H=1,T=2,order=1,stayer=1,occasion=FALSE),
  data.frame(K=2,H=1,T=3,order=1,stayer=1,occasion=FALSE),
  data.frame(K=3,H=1,T=3,order=1,stayer=1,occasion=FALSE),
  data.frame(K=2,H=1,T=3,order=1,stayer=1,occasion=TRUE),
  data.frame(K=3,H=1,T=3,order=1,stayer=1,occasion=TRUE),
  data.frame(K=2,H=1,T=4,order=1,stayer=1,occasion=TRUE),
  data.frame(K=2,H=2,T=3,order=1,stayer=0,occasion=TRUE),
  data.frame(K=2,H=2,T=4,order=1,stayer=0,occasion=TRUE),
  data.frame(K=3,H=2,T=3,order=1,stayer=0,occasion=TRUE),
  data.frame(K=2,H=1,T=3,order=2,stayer=0,occasion=FALSE),
  data.frame(K=3,H=1,T=3,order=2,stayer=0,occasion=FALSE),
  data.frame(K=2,H=2,T=3,order=2,stayer=0,occasion=FALSE),
  data.frame(K=2,H=2,T=4,order=2,stayer=0,occasion=FALSE),
  data.frame(K=2,H=1,T=3,order=2,stayer=1,occasion=FALSE),
  data.frame(K=2,H=1,T=4,order=2,stayer=1,occasion=FALSE),
  data.frame(K=2,H=1,T=4,order=2,stayer=0,occasion=TRUE))
grid$counted <- mapply(count, grid$K, grid$H, grid$T, grid$order, grid$stayer, grid$occasion)
grid$saturated <- grid$K^grid$T - 1
grid$rank <- mapply(rank_at, grid$K, grid$H, grid$T, grid$order, grid$stayer, grid$occasion)
grid$excess <- grid$counted - grid$rank
print(grid)
