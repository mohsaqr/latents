# multilpa 0.11.0

This release is a correctness release. An independent review of 0.10.0 found
that several advertised paths returned plausible numbers that were wrong, rather
than failing. Every one of them is fixed below, each with a regression test that
reproduces the original defect. Nothing here adds a model family.

## Silently wrong results, now correct

### A covariate's units no longer change the fit

`fit_covariates()` ran its membership M-step on the coefficients as supplied.
BFGS stops on a *relative* tolerance measured in the coordinates it is given, so
a predictor in small units stalled the search at a point where the score was
still large, and `optim()` reported success. Rescaling one covariate by `1e-7`
moved the log likelihood from `-620.1281` to `-730.5748` and drove the slope to
zero, with no warning. The search now runs on root-mean-square-scaled
coordinates and convergence is judged by a unit-free score, not by the
optimiser's return code. Estimates are returned on your scale, unchanged:

```r
# identical to machine precision at every scale
fit_covariates(data, "y", "g", 2, profile_covariates = "z")
```

### Staged categorical measurement attaches to the right item

A held measurement was passed between stages as an unlabelled, positional list,
while the compatibility check compared only the *sorted set* of indicator names
-- which a permutation passes. Reordering `categorical = c("a", "b")` to
`c("b", "a")` swapped whole response distributions (`-587.5909` to `-671.5937`);
reversing a factor's levels swapped categories within an item (`-1047.9689`).
`starting_values()` now keeps the indicator names and category labels, and
`multilpa()` aligns on them. An encoding that cannot be aligned raises
`multilpa_bad_start` or `multilpa_bad_stage` instead of returning a number.
An unlabelled start is still read positionally, so stored starts keep working.

### Bootstrap comparisons keep the constraints they were given

`bootstrap_lrt()` accepted fixed-measurement models but dropped `fixed` and the
held values from every refit, so the observed statistic compared constrained
fits while the null distribution came from unconstrained ones -- a 1-versus-3
parameter statistic referred to a 9-versus-11 parameter null. Refits now carry
the constraint, and a pair whose constrained nesting cannot be established is
refused with `multilpa_bad_nesting`.

### Transition counts no longer overflow

Expected transition counts exponentiated a scaled forward factor and an unscaled
backward factor separately, so on a long sequence one overflowed to `Inf` where
the other underflowed to `0`. A converged 600-occasion fit reported `NaN` for
every count, and fitting 1,800-occasion sequences failed outright. The complete
pair log probability is now formed before exponentiating, which cannot overflow.
Counts agree with explicit path enumeration to `1.8e-15`.

### `descriptives()` no longer invents or replaces values

An `NA` stratifier produced `NA` row selectors, and `x[NA, ]` fabricates a row
rather than dropping one: every stratum gained a phantom missing observation
while the real row vanished. Missing strata are now a labelled stratum, and the
per-stratum counts are asserted to account for every input row. Separately, the
profile assignment was written into the same frame as the indicators, so an
indicator named `profile` or `group_class` was overwritten by its own class
labels. The assignment is now kept beside the frame, never in it.

### Cluster-robust standard errors are refused when they cannot exist

`three_step()` and `r3step()` summed estimating-equation contributions within
groups without requiring enough independent groups. Those contributions sum to
zero at the estimate -- that sum *is* the stationarity condition -- so a single
group gave a standard error of `5.9e-17` and a confidence interval of zero
width, and `r3step()`'s robust path reported p-values of exactly zero. Both now
require more independent groups than the quantities reported, raising
`multilpa_too_few_groups`. `three_step(vcov_type = "independent")` offers an
explicit, labelled unclustered alternative.

### `diagnostics()` no longer ignores the argument you passed

`diagnostics()` and `report()` documented that `...` reached the verbs they
gather, but discarded it. Because `by` is a real argument of
`bivariate_residuals()`, `diagnostics(fit, by = "overall")` silently returned
per-profile residuals under the name you asked to pool. Both verbs now take `by`
explicitly and refuse any argument they cannot forward.

### A bootstrap comparison stops withholding its p-value

`bootstrap_lrt()` judged a replicate invalid if its statistic fell below a fixed
`-1e-5`. But the statistic is `2 * (ll_alt - ll_null)` on log likelihoods of
order 10^3 to 10^4, while EM stops on a *relative* tolerance, and under the null
the alternative converges to the null solution -- so every replicate sits at the
boundary with noise of about `2 * tol * |log likelihood|`. On
`school_engagement` that is 7.6e-05, seven times the window. Twelve replicates
that all converged had four of them reported as "Nonconvergence or reversed
likelihood" and the p-value withheld. The window now scales with `tol` and the
likelihood, so the same comparison returns 12 of 12 valid and a p-value.

### An ordinary converged fit stops being told it has not converged

`parameter_inference()` cut a fixed `0.01` on the scaled score. On 720
observations at the default `tol = 1e-8` that fires on a fit which converged
with all 20 starts replicating. The criterion is now the displacement the score
implies, expressed in the estimate's own standard errors: a score `g` against
information `I` moves the estimate by `g / I`, and the standard error is
`sqrt(1 / I)`, so `g * SE` is that displacement and is dimensionless. The
default-tolerance fit sits 0.17% of a standard error from stationarity and is
now quiet; a fit at `tol = 1e-4` sits 5% away and still warns, now saying so in
those words.

### Misaligned data is refused, not scored

`assignments()` documented that it checked row alignment but compared only the
row count, and `bivariate_residuals()` trusted row order. A reversed frame was
paired with the wrong posteriors and returned as a result. Both now compare
every column the fit recognises and raise `multilpa_bad_inference_data` on a
mismatch; a frame sharing no column with the fit warns
`multilpa_unverified_alignment` rather than passing in silence.

## Now supported

* `parameter_inference()`, `vcov()` and `confint()` work on fixed and staged
  fits, reporting standard errors conditional on the held measurement -- which
  `?fit_staged` already promised. They previously failed with
  `length(theta) == object$n_parameters is not TRUE`. Held blocks enter the
  likelihood at their held values and contribute no row to the information
  matrix; the conditional likelihood matches an independent enumeration to
  `5.7e-14` and the standard errors match an independent numerical information
  matrix to `1.3e-07` relative.
* `assignments()` and `descriptives(by = "profile")` work on a random-intercept
  fit. Both were documented to, but the fit stored no modal assignment, so the
  first failed with an unclassed `cbind()` error.
* `multilpa_plot_types()` lists `"random_intercepts"`, which
  `plot.multilpa_random_intercept()` has always accepted.
* `fit_transitions(max_iter = 0)` returns a usable fit instead of failing with
  "attempt to set an attribute on NULL".
* `plot()` draws `what = "entropy"` and `what = "posteriors"` for covariate and
  random-intercept fits, so `plot(diagnostics(fit))` works for every family.

## Behaviour changes

* `max_iter = 0` with a supplied `start` now evaluates that start alone and
  ignores `n_starts`, which is what `?multilpa` promised. It previously scored
  all random starts and returned the best of them.
* Restart selection now breaks ties deterministically rather than by
  `which.max()`, which was choosing between equally valid optima on
  floating-point noise.
* A separated covariate fit is now reported as unconverged rather than
  converged, so `parameter_inference()` refuses it instead of returning Wald
  intervals around an unbounded coefficient.
* Every warning the package raises now carries a documented condition class, so
  a simulation loop can muffle the qualification it expects and let the rest
  through. New classes: `multilpa_extreme_coefficients`,
  `multilpa_quadrature_check`, `multilpa_empty_transition_row`,
  `multilpa_unverified_alignment`, `multilpa_no_free_parameters`,
  `multilpa_held_parameter`.
* `three_step()` gains `vcov_type`; `diagnostics()` and `report()` gain `by`;
  `plot.multilpa_covariates()` and `plot.multilpa_random_intercept()` gain
  `"entropy"` and `"posteriors"`. All are additive and keep positional calls.
* `multilpa_bootstrap_lrt` objects carry `fixed`, and
  `as.data.frame(what = "test")` gains a `fixed` column.
* `parameter_inference()` on a fixed fit reports no row for a held block and
  carries a `fixed` attribute. The table is on the natural scale, where every
  probability in a set is reported, so it has `n_parameters` rows plus one for
  each reference category the estimation scale drops. `vcov()` is
  `n_parameters` square on the unconstrained scale; on the natural scale it is
  correspondingly larger and singular.
* `descriptives(by = )` can return one additional row: the labelled missing
  stratum.
* `profile_prevalence` for a transition fit is read from the final expectation
  rather than the last M-step. These agree at convergence; the new value is the
  correct one for an unconverged or zero-iteration fit.
* `starting_values()` keeps names and column labels on categorical response
  blocks, and `as.data.frame(what = "responses")` now reports those labels
  rather than positions, matching what the fitted object's own response table
  reports for the same quantity. A start carrying no labels is still read
  positionally and its tidy view still reports positions.
* `assignments()` raises `multilpa_bad_inference_data` rather than
  `multilpa_bad_nesting` for a row-count mismatch, which is a data contract
  failure, not a model-nesting one.
* Thirteen further user-reachable guards in `bootstrap_lrt()` and
  `fit_transitions()` now raise their catalogued class instead of a bare
  `simpleError`, so `?"multilpa-conditions"` is true of them.
* A `multilpa(fixed = )` fit records `n_parameters_with_measurement`, which only
  `fit_staged()` set before. `print()` consequently reported the same count
  twice on such a fit; it now reports the two counts it names.
* `parameter_inference()` results carry `score_displacement` beside
  `scaled_score`, the quantity the convergence warning is judged on.

## Infrastructure

* GitHub Actions check the package on macOS, Windows and Linux across release,
  devel and oldrel, plus an explicit R 4.1 job for the declared `Depends` floor;
  a second workflow runs `lintr`.
* `DESCRIPTION` gains `URL` and `BugReports`.

# multilpa 0.10.0

## New features

### A fit carries the data it was built from

Post-fit verbs no longer ask for data the fit is already holding. The model
frame is rebuilt from the uncentred indicators, the categorical codes and their
levels, the identifier and the occasion — all of which were stored already — so
this costs no extra memory, and the round trip is exact for every model family.

```r
parameter_inference(fit)          # was parameter_inference(fit, students)
bivariate_residuals(fit)
vcov(fit); confint(fit)
bootstrap_lrt(smaller, larger)
as.data.frame(fit, what = "data") # the columns the model was fitted to
```

Supplying `data` still works and is still the stronger check, since a frame that
does not reproduce the fitted likelihood is rejected as before.
`three_step()` and `r3step()` keep requiring `data`, because an outcome or a
covariate is a column the model never saw. The `multilpa_data_required`
condition is retired: there is now something sensible to default to.

### Three convenience verbs

* `descriptives(x, vars, id)` — one row per variable with n, missing, mean, sd,
  range, and the **intraclass correlation**: the share of variance lying between
  groups. An ICC near zero says the groups do not differ, so a model built to
  tell them apart has nothing to find. It accepts a data frame before fitting or
  a fit afterwards, and `by = "profile"` splits it by assigned profile.
* `diagnostics(x)` — entropy, classification quality, average posteriors and
  bivariate residuals in one call, as a classed object with `print()`,
  `as.data.frame(x, what = )` and `plot()`. `plots = TRUE` draws them too. A
  model family that cannot supply a table says so rather than failing.
* `report(x)` — summary, descriptives, diagnostics and every plot the fit
  supports, in one call. It skips views a given fit has no ingredients for
  rather than stopping at the first refusal.

### `assignments()`

Every observation with the profile it was assigned to, the group class, and its
posteriors -- optionally with a frame of your own columns kept alongside:

```r
labelled <- assignments(fit, data = school_engagement)
xtabs(~ profile + engaged, data = labelled)
```

Comparing an assignment with anything else means putting them in the same row,
and doing that with two separate objects assumes they share an order. This verb
owns the alignment: a frame with the wrong number of rows is refused, and a
column the assignments would overwrite is an error rather than a replacement.

It is a verb rather than another `as.data.frame(what = )` value because it
computes a join with a frame the caller supplies. The line: `as.data.frame()`
represents what a fit already holds, and anything bringing in columns the model
never saw gets its own verb with its own documented arguments.

### `plot(fit, what = "bars")` draws its intervals without being handed data

The fit carries the columns it was built from, so the 95% intervals are drawn by
default; `data` is now only an override. A model family whose standard errors
are not implemented gets bars without whiskers rather than an error.

`report()` likewise asks each `plot()` method what views it accepts instead of
assuming the full catalogue, and skips a view the method refuses -- naming what
it could not draw rather than dropping it silently.

## Bug fixes

* The intraclass correlation is now `1` rather than `NA` when a variable has no
  within-group variance. That is a well-defined case — every observation equals
  its group's own value — not a degenerate one.

# multilpa 0.9.0

## Breaking changes

### `information_criteria()` returns the shape people report

The default is now one row, with one column per criterion — the shape a
model-comparison table is published in, and the same column names
`enumerate_classes()` uses, so a chosen model's row sits under the grid it was
chosen from without renaming anything.

```r
information_criteria(fit)
#>   log_likelihood n_parameters      aic      kic bic_groups bic_individual ...
#> 1      -3778.259           15 7586.518 7604.518   7617.933       7655.207 ...

information_criteria(fit, format = "long")   # the previous shape
```

* The `penalty` column has been removed. No convention reports it — not Mplus,
  not `tidyLPA::get_fit()`, not mclust — and it is the difference between two
  numbers the table already carries.
* `definitions = TRUE` describes one criterion per row, so it now requires
  `format = "long"` and raises `multilpa_bad_argument` otherwise.
* `as.data.frame(fit, what = "information_criteria")` follows the verb, and its
  `format` argument now actually reaches it. Previously `format` was captured by
  the accessor's own signature and never forwarded.
* `convention` remains `NA_character_` for `deviance`, `aic` and `kic`, now
  documented as meaning *the question does not arise* rather than *a value is
  missing*.

## New features

* Two example datasets ship with the package, so the first example runs as
  written: `school_engagement` (720 students in 60 schools) and
  `engagement_panel` (120 students over four waves). Both are simulated and both
  carry the kind each row was generated from — `engaged` and `state` — so a
  fitted model can be checked against what produced it.
* `vignette("multilpa")` is rebuilt on those datasets and gains a section
  cross-tabulating the fitted profiles against the truth column.

# multilpa 0.8.0

## Breaking changes

### Argument names now follow the sibling package `Nestimate`

No deprecation shim: the old spellings are gone and raise
`unused argument`.

| was | is | where |
|---|---|---|
| `object` | `x` | first formal of this package's own verbs |
| `indicators` | `vars` | `multilpa()`, `fit_covariates()`, `fit_staged()`, `fit_transitions()`, `fit_random_intercept()`, `enumerate_classes()` |
| `group` | `id` | the same six |
| `level_ci` | `ci_level` | `three_step()`, `r3step()` |
| `p_adjust` | `adjust` | `parameter_inference()`, `three_step()`, `r3step()` |
| `n_boot` | `iter` | `bootstrap_lrt()` |
| `profiles`, `group_classes` | `n_profiles`, `n_group_classes` | `enumerate_classes()` |

```r
# before
multilpa(students, indicators = vars, group = "school", n_profiles = 2)
# after
multilpa(students, vars = vars, id = "school", n_profiles = 2)
```

The S3 methods for the base generics `summary()`, `coef()`, `vcov()`,
`confint()`, `logLik()` and `nobs()` keep `object`, because those generics own
that name.

**Tidy output columns did not change.** Every posterior, sequence and
classification table keeps its `group` column, `level` keeps the value
`"group"`, and fitted objects keep `group_classes` and `n_group_classes`.

## New features

* `multilpa_plot_types()` lists every view the `plot()` methods accept, with the
  group each belongs to and what it answers. It carried `@export` in 0.7.0 but
  was missing from `NAMESPACE`, so it could not be called.
* `plot()` gains four views: `"bars"` (profile means as grouped bars, with 95%
  intervals when `data` is supplied), `"heatmap"` (means as deviations from each
  indicator's grand mean), `"entropy"` and `"posteriors"` (per-case
  classification quality as ridges).

## Bug fixes

* A refused plot no longer leaves the graphics device unusable. Restoring `mfg`
  as part of `par(no.readonly = TRUE)` switches `new` on as a documented side
  effect, so a view that validated and raised before drawing left every
  subsequent plot ready to overlay onto a blank panel.
* Direct labels are no longer clipped at the right edge on narrow panels. The
  gap was measured in user units while the margin reserved for it was measured
  in inches.
* Point area in `plot(fit)` encodes prevalence absolutely rather than being
  min-max normalised within the plot, which had mapped any spread onto the full
  size range and so carried no information.

# multilpa 0.7.0

* Every public verb returns a tidy `data.frame`. See `CHANGES.md` in the source
  repository for the full sweep.
