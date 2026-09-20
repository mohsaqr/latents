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
