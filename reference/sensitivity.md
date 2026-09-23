# Does the solution survive a different seed?

A latent profile model is fitted by EM from random starts, so what comes
back is a local maximum reached from one particular set of them.
`n_starts` reports how many starts *within a single seed's stream*
reached the best likelihood; it cannot say whether a different stream
would have found a different mode. This verb refits the model under each
of several seeds and reports what changed.

## Usage

``` r
sensitivity(
  x,
  data = NULL,
  seeds = 1:10,
  n_starts = NULL,
  max_iter = NULL,
  tol = NULL,
  tolerance = NULL
)
```

## Arguments

- x:

  A fitted `multilpa` model to use as the reference. A
  [`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
  result varies the second-stage starts while holding its first-stage
  measurement solution fixed. A directly fitted model with measurement
  blocks held is refused: its original free starting values are not
  retained, so replaying the same seed would not replay the same search.

- data:

  Optional. The data the model was fitted to, in its original row order.
  A fit carries the columns it was built from, so this is only needed to
  override them. Shared identifier and indicator columns are checked
  against the fit before assignment agreement is calculated.

- seeds:

  Seeds to refit under: at least two distinct whole numbers of either
  sign, as [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts.
  Anything else raises `latents_bad_argument`. The reference fit's own
  seed may be among them, in which case that row reproduces it and is
  the arithmetic check that the refit is the same model.

- n_starts:

  Starts per refit. Defaults to the number the reference fit used, so
  each seed gets the same search effort it did.

- max_iter, tol:

  Passed to each refit; default to the reference fit's.

- tolerance:

  Two log likelihoods within this of each other count as the same
  optimum. Defaults to the reference fit's own replication tolerance,
  which is what it used to decide whether its starts agreed.

## Value

A base `data.frame`, one row per seed, with columns `seed`,
`log_likelihood`, `converged`, `iterations`, `optimum` (1 for the best
maximum found across the seeds, 2 for the next distinct one, and so on),
`best` (whether this seed reached optimum 1) and `agreement` (the
proportion of observations given the same profile as the reference fit,
after label alignment). A failed refit has `NA` throughout except
`seed`.

## Details

Three things are worth reading off the result. Whether every seed
reached the same maximised log likelihood, which is what `optimum`
counts. Whether the seeds that reached it converged. And how many
observations were assigned to a different profile, which is what
`agreement` measures, after the arbitrary profile labels of each refit
have been matched to the reference fit's – two fits of the same mixture
can be identical and still number their profiles differently, so
comparing labels directly would report disagreement that is not there.

A seed whose refit fails contributes a row with `NA` estimates rather
than being dropped silently, and a `latents_sensitivity_dropped` warning
names how many failed, so the table is never quietly shorter than
`seeds`.

## References

Hipp, J. R., & Bauer, D. J. (2006). Local solutions in the estimation of
growth mixture models. *Psychological Methods, 11*, 36–53.

## See also

[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
to vary the number of classes rather than the seed, and
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
for uncertainty about the parameters of a single solution.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = c(rnorm(60, -1.5), rnorm(60, 1.5)),
  score_b = c(rnorm(60, -1), rnorm(60, 1))
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 3, seed = 1)
sensitivity(fit, seeds = 1:3, n_starts = 3)
#>   seed log_likelihood converged iterations optimum best agreement
#> 1    1      -396.1633      TRUE         13       1 TRUE         1
#> 2    2      -396.1633      TRUE         13       1 TRUE         1
#> 3    3      -396.1633      TRUE         13       1 TRUE         1
```
