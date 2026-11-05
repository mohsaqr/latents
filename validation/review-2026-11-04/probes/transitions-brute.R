# Independent brute-force check of the likelihood, posteriors and expected
# counts: enumerate every latent path of every person in plain base R, reading
# parameters off the fit and building covariate designs by hand from the data.
suppressMessages(devtools::load_all("/Users/mohammedsaqr/Documents/Github/latents", quiet = TRUE)); pdf(NULL); set.seed(1)

softmax <- function(eta) { eta <- eta - max(eta); exp(eta) / sum(exp(eta)) }

# Person-level brute force. `grid_of` gives the slot of each row.
# init(i, h) -> K-vector; move(i, h, s, prev, cur, to_slot) -> K-vector for next state.
brute <- function(d, vars, K, H, shares, emit, init, move, slots, order = 1L, weights = NULL) {
  ids <- unique(d$id)
  out <- lapply(seq_along(ids), function(i) {
    rows <- which(d$id == ids[i])
    s <- slots[rows]; L <- max(s)
    paths <- as.matrix(expand.grid(rep(list(seq_len(K)), L)))
    logp <- vapply(seq_len(H), function(h) {
      apply(paths, 1L, function(path) {
        lp <- log(init(i, h)[path[1]])
        if (L >= 2L) for (t in 2:L) { # small scratch loop in a probe
          lp <- lp + log(move(i, h, t, if (t >= 3L) path[t-2] else NA, path[t-1], path[t])[1])
        }
        for (r in rows) lp <- lp + emit(r, path[slots[r]])
        lp
      })
    }, numeric(nrow(paths)))
    logp <- matrix(logp, nrow(paths), H)
    joint <- sweep(logp, 2L, log(shares), "+")
    m <- max(joint); ll <- m + log(sum(exp(joint - m)))
    post <- exp(joint - ll)
    # posterior per row by profile
    row_post <- t(vapply(rows, function(r) vapply(seq_len(K), function(k)
      sum(post[paths[, slots[r]] == k, ]), 0), numeric(K)))
    counts <- array(0, c(K, K, H))
    if (L >= 2L) for (t in 2:L) for (a in seq_len(K)) for (b in seq_len(K)) for (h in seq_len(H))
      counts[a, b, h] <- counts[a, b, h] + sum(post[paths[, t-1] == a & paths[, t] == b, h])
    list(ll = ll, rows = rows, row_post = row_post, counts = counts,
         class_post = colSums(post), init_post = vapply(seq_len(H), function(h)
           vapply(seq_len(K), function(k) sum(post[paths[,1]==k, h]), 0), numeric(K)))
  })
  rp <- matrix(0, nrow(d), K)
  for (o in out) rp[o$rows, ] <- o$row_post
  list(ll = sum(vapply(out, `[[`, 0, "ll")), row_post = rp,
       counts = Reduce(`+`, lapply(out, `[[`, "counts")),
       class_post = matrix(t(vapply(out, `[[`, numeric(H), "class_post")), length(out), H),
       init_counts = Reduce(`+`, lapply(out, `[[`, "init_post")))
}

simulate <- function(n = 40, T = 4, K = 2, seed = 1) {
  set.seed(seed)
  P <- matrix(c(.8,.2,.3,.7), 2, byrow = TRUE)
  rows <- lapply(seq_len(n), function(i) {
    s <- integer(T); s[1] <- sample(2, 1, prob = c(.6,.4))
    for (t in 2:T) s[t] <- sample(2, 1, prob = P[s[t-1], ])
    data.frame(id = i, time = seq_len(T), y1 = rnorm(T, c(0, 2)[s]), y2 = rnorm(T, c(1, -1)[s], .8))
  })
  d <- do.call(rbind, rows)
  d$z <- rnorm(nrow(d)); d$w <- rep(rnorm(n), each = T)
  d
}

report <- function(label, bf, ll, post, counts = NULL) {
  cat(sprintf("%-46s dLL %.2e  max|dpost| %.2e  max|dcount| %s\n", label, ll - bf$ll,
              max(abs(unname(post) - bf$row_post)),
              if (is.null(counts)) "NA" else sprintf("%.2e", max(abs(unname(counts) - bf$counts)))))
}

base <- simulate()
# Gaps: person-specific skipped positions, plus dropout, shuffled rows
g <- base
drop <- with(g, (id %% 3 == 0 & time == 2) | (id %% 5 == 0 & time == 3) | (id %% 7 == 0 & time == 4))
g <- g[!drop, ]
late <- g$id %% 11 == 0 & g$time == 1   # late entry
g <- g[!late, ]
g <- g[sample(nrow(g)), ]
stopifnot(any(!drop))
vars <- c("y1", "y2")

slots_grid <- function(d) { pos <- sort(unique(d$time)); p <- match(d$time, pos); ave(p, d$id, FUN = function(v) v - min(v) + 1L) }
slots_obs  <- function(d) ave(d$time, d$id, FUN = rank)

emit_of <- function(d, means, vars_mat) function(r, k) sum(dnorm(unlist(d[r, vars]), means[k, ], sqrt(vars_mat[k, ]), log = TRUE))

# 1 homogeneous engine on a grid with gaps, one class and two classes
for (H in 1:2) for (occ in c("grid", "observed")) {
  f <- suppressWarnings(lta(g, vars, "id", 2, "time", n_group_classes = H, n_starts = 2, seed = 1, occasions = occ, max_iter = 30))
  s <- if (occ == "grid") slots_grid(g) else slots_obs(g)
  bf <- brute(g, vars, 2, H, f$group_probabilities, emit_of(g, f$means, f$variances),
              function(i, h) f$initial_probabilities[h, ],
              function(i, h, t, p, a, b) f$transition_probabilities[a, b, h], s)
  report(sprintf("homogeneous %s H=%d", occ, H), bf, f$log_likelihood, f$subject_posteriors, f$transition_counts)
  cat(sprintf("   class post max diff %.2e; prevalence vs row post %.2e\n",
     max(abs(unname(f$group_posteriors) - bf$class_post)),
     if (H == 1) max(abs(colMeans(bf$row_post) - f$profile_prevalence[1, ])) else NA))
}

# 2 general engine: occasion transitions with gaps on grid
f <- suppressWarnings(lta(g, vars, "id", 2, "time", transitions = "occasion", occasions = "grid", n_starts = 2, seed = 1, max_iter = 30))
tc <- f$transition_coefficients
s <- slots_grid(g)
bf <- brute(g, vars, 2, 1, f$group_probabilities, emit_of(g, f$measurement[[1]]$means, f$measurement[[1]]$variances),
  function(i, h) softmax(c(f$initial_coefficients[1,1,h], 0)),
  function(i, h, t, p, a, b) { other <- setdiff(1:2, a); e <- numeric(2); e[other] <- tc[t-1, 1, a, h]; softmax(e)[b] }, s)
report("general occasion grid gaps", bf, f$log_likelihood, f$subject_posteriors)
print(dimnames(tc)[[1]])

# 3 general engine, second order, grid gaps, 2 classes
f <- suppressWarnings(lta(g, vars, "id", 2, "time", order = 2, occasions = "grid", n_group_classes = 2, n_starts = 2, seed = 1, max_iter = 30))
tc <- f$transition_coefficients; t2 <- f$second_order_coefficients
bf <- brute(g, vars, 2, 2, f$group_probabilities, emit_of(g, f$measurement[[1]]$means, f$measurement[[1]]$variances),
  function(i, h) softmax(c(f$initial_coefficients[1,1,h], 0)),
  function(i, h, t, p, a, b) {
    other <- setdiff(1:2, a); e <- numeric(2)
    e[other] <- if (t == 2L) tc[1, 1, a, h] else t2[1, 1, (p - 1L) * 2L + a, h]
    softmax(e)[b] }, s)
report("general order 2 grid gaps H=2", bf, f$log_likelihood, f$subject_posteriors)

# 4 stayer, grid gaps, occasion measurement, initial covariate
f <- suppressWarnings(lta(g, vars, "id", 2, "time", mover_stayer = TRUE, occasions = "grid", measurement = "occasion", initial_covariates = "w", n_starts = 2, seed = 1, max_iter = 30))
tc <- f$transition_coefficients
first_rows <- vapply(unique(g$id), function(i) { r <- which(g$id == i); r[which.min(g$time[r])] }, 1L)
occ_slot <- s
emit_occ <- function(r, k) { m <- f$measurement[[occ_slot[r]]]; sum(dnorm(unlist(g[r, vars]), m$means[k, ], sqrt(m$variances[k, ]), log = TRUE)) }
bf <- brute(g, vars, 2, 2, f$group_probabilities, emit_occ,
  function(i, h) { b <- f$initial_coefficients[, 1, h]; softmax(c(b[1] + b[2] * g$w[first_rows[i]], 0)) },
  function(i, h, t, p, a, b) { if (h == 2L) return(as.numeric(a == b)); other <- setdiff(1:2, a); e <- numeric(2); e[other] <- tc[1, 1, a, h]; softmax(e)[b] }, s)
report("stayer + occ measurement + init cov, grid gaps", bf, f$log_likelihood, f$subject_posteriors)

# 5 transition covariate observed layout, destination reading, factor reference coding, K = 3
d3 <- simulate(n = 25, T = 3)
d3$f <- factor(c("b", "a", "c")[(d3$id %% 3) + 1L], levels = c("c", "a", "b"))
d3 <- d3[sample(nrow(d3)), ]
f <- suppressWarnings(lta(d3, vars, "id", 3, "time", transition_covariates = c("z", "f"), n_starts = 2, seed = 1, max_iter = 30))
print(dimnames(f$transition_coefficients)[[1]])
tc <- f$transition_coefficients
s3 <- slots_obs(d3)
row_at <- function(i, t) which(d3$id == unique(d3$id)[i] & d3$time == t)
bf <- brute(d3, vars, 3, 1, 1, emit_of(d3, f$measurement[[1]]$means, f$measurement[[1]]$variances),
  function(i, h) softmax(c(f$initial_coefficients[1, , h], 0)),
  function(i, h, t, p, a, b) {
    r <- row_at(i, t)
    z <- c(1, d3$z[r], d3$f[r] == "a", d3$f[r] == "b")
    e <- numeric(3); others <- setdiff(1:3, a)
    e[others] <- as.vector(t(tc[, , a, h]) %*% z); softmax(e)[b] }, s3)
report("K=3 transition covariates, factor levels c/a/b", bf, f$log_likelihood, f$subject_posteriors)
cat("n params check K=3 cov:", f$n_parameters, " hand: means 3*2 + var 3*2 + init 2 + trans 4 columns*2 dest*3 origins =", 6+6+2+4*2*3, "\n")
