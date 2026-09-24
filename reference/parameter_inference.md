# Tidy inference for a covariate fit

Standard errors, tests and intervals for every free parameter of
`multilpa(profile_covariates = )`, at all three levels at once: the
measurement model, the profile logits and the group-class logits.
Without them a membership coefficient cannot be reported, because
nothing distinguishes a real effect from separation.

Computes ordinary maximum-likelihood standard errors from the inverse
observed Hessian and delta-method Wald confidence intervals. These are
local, asymptotic estimates conditional on the selected numbers of
classes. They are not robust sandwich errors and do not account for
model selection. Mixture boundary solutions and unidentified Hessians do
not admit this calculation.

## Usage

``` r
# S3 method for class 'multilpa_covariates'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL
)

parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL
)

# S3 method for class 'multilpa'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL
)
```

## Arguments

- x:

  A converged `multilpa` fit with inactive variance bounds.

- data:

  Optional. The data frame the model was fitted to; when omitted it is
  rebuilt from the indicators, identifiers and occasions the fit stores,
  which round-trip exactly. Supplying it is the stronger check that the
  caller still holds that frame. A supplied frame is checked exactly
  against the stored indicators and identifiers; for an older fit
  without stored indicators, the likelihood provides a weaker
  consistency check.

- level:

  Confidence level strictly between zero and one.

- step:

  Positive finite-difference step for the analytic score derivative.

- vcov_type:

  `"observed"` inverts the observed information. `"robust"` returns the
  Huber-White sandwich `A^-1 B A^-1`, where `B` accumulates the outer
  product of the per-group scores. Groups are the independent units, so
  robust errors relax the assumption that the Gaussian within-group
  model is correctly specified. They do not relax the assumption that
  groups are independent.

- adjust:

  Multiplicity correction applied to `p_value` to produce `p_adjusted`,
  one of the methods
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) accepts.
  The default `"none"` leaves the two columns equal: a correction
  changes what a p-value means, so it is applied only when it is asked
  for, and the method that was applied is recorded in the `adjust`
  attribute. The correction is taken over the tests the table actually
  reports; bounded parameters carry no test and do not count towards the
  family.

- method:

  How uncertainty is measured. `"wald"`, the default, inverts the
  observed information in the estimation coordinates, and is available
  for the four covariance structures those coordinates can express: EEI,
  VVI, EEE and VVV. `"bootstrap"` resamples *groups* with replacement,
  refits inside the same covariance family, undoes the relabelling each
  refit comes back with, and reports the percentile interval and the
  standard deviation of the replicates. It is available for all fourteen
  structures, because it needs no coordinate chart for the constraint.
  Groups are the resampling unit rather than rows, so the interval
  carries the same independence assumption as `vcov_type = "robust"` and
  not the stronger one `"observed"` makes.

- iter:

  Number of resamples when `method = "bootstrap"`, at least two. The
  default 199 is enough for a standard error; a 95% percentile interval
  is steadier at 999 and above.

- n_starts, max_iter, tol:

  Passed to each bootstrap refit. Each replicate is a fresh fit, so
  `n_starts` buys the same protection against a local maximum here as it
  does in
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  at the same multiple of the cost.

- seed:

  Optional seed for the resampling: any whole number
  [`set.seed()`](https://rdrr.io/r/base/Random.html) accepts; a fraction
  raises `latents_bad_argument`. The stream is restored on exit, so a
  seeded call leaves the caller's random state exactly as it found it.

## Value

A base `data.frame` with one row per free parameter and the same
columns, in the same order, that `parameter_inference()` returns for
every other fitted class: `level` (`"measurement"`, `"profile"` or
`"group"`), `outcome` (which profile or group class the coefficient
predicts, or which profile a measurement parameter belongs to), `term`,
`parameter` (`"mean"`, `"variance"`, `"covariance"` or `"coefficient"`),
`estimate`, `standard_error`, `statistic`, `p_value`, `p_adjusted`,
`conf_low` and `conf_high`. Variances and covariances are reported in
their natural units, with standard errors carried through the delta
method from the log and log-Cholesky coordinates on which they are
estimated.

A base `data.frame` with one row per reported parameter and the columns
`level` (`"measurement"`, `"profile"` or `"group"`), `outcome`, `term`,
`parameter`, `estimate`, `standard_error`, `statistic`, `p_value`,
`p_adjusted`, `conf_low` and `conf_high`. Estimates are on the natural
scale: variances in their own units, probabilities as probabilities. A
parameter whose null value is on the boundary of its space – a variance,
a probability, a residual variance on the covariance diagonal – carries
`NA` for `statistic`, `p_value` and `p_adjusted`, because a Wald test
against that boundary is not a question worth asking; the interval still
is. All class probabilities are reported; their sum constraints make the
natural covariance singular. Wald intervals are not clipped to parameter
bounds.

A fit made with `fixed` reports only the parameters it estimated: a held
measurement block is a constant of this likelihood, so it contributes no
row here and no row or column to the information matrix. The table is on
the natural scale, where every class probability is reported, so it has
one row per estimated natural coefficient: `n_parameters` rows plus one
for each set of probabilities whose reference category the estimation
scale drops. The held values still enter the likelihood at the values
they were held at, so every standard error, interval and p-value is
conditional on that measurement solution and does not propagate its
uncertainty. See
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
for what that conditioning means.

With `method = "bootstrap"` the table has the same columns and the same
rows, `estimate` is still the fitted value, `standard_error` is the
standard deviation of the replicates and `conf_low`/`conf_high` are the
percentile interval. `statistic`, `p_value` and `p_adjusted` are `NA`
throughout: putting a normal approximation back on top of the replicates
the interval was read from is the assumption this path exists to avoid,
so the interval is the inference. Its attributes are `covariance` (the
covariance of the replicates, which
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns), `level`,
`method`, `iter`, `n_valid`, `replicates` (the kept replicates, one row
each), `messages` (why a resample was dropped, `NA` where it was not),
`structure` and `fixed`; `covariance_unconstrained` is `NULL`, because
the estimation-scale chart is what this path does without.

Diagnostics of the fit as a whole travel as attributes rather than as
columns repeated down every row: `covariance` and
`covariance_unconstrained` (natural and estimation-scale covariance of
the estimates), `hessian`, `gradient`, `scaled_score`,
`condition_ratio`, `level`, `step`, `vcov_type`, `fixed` (the
measurement blocks this fit held,
[`character()`](https://rdrr.io/r/base/character.html) when none) and
`adjust`. When `vcov_type` is `"robust"`, `group_scores` holds the
groups-by-parameters score matrix and `scaling_correction` the MLR
scaling correction factor `tr(A^-1 B) / q`, used for scaled
likelihood-ratio difference tests.

## Details

Class-membership coefficients are on the multinomial logit scale with
the final profile and the final group class as references, so a
coefficient is a log odds against that reference. The tests are Wald
tests and inherit the usual caveat: they are unreliable for a
coefficient driven to the boundary by separation, which the size of the
estimate and its standard error together will reveal.

## Conditions

`latents_no_converge` for an unconverged fit, `latents_boundary_fit`
when a variance, a response probability or a mixing probability sits at
its bound, `latents_singular_information` when the observed information
cannot be inverted, `latents_no_free_parameters` when `fixed` held every
parameter, `latents_bad_inference_data` when a supplied `data` does not
reproduce the fit, and `latents_too_few_groups` for
`vcov_type = "robust"` with fewer independent groups than reported
parameters. A fit whose score is still far from zero is reported with a
`latents_unconverged` warning rather than refused.

`method = "bootstrap"` adds `latents_unsupported_inference` for a fit
that holds a measurement block — the held values came from another fit
and resampling these data does not resample them — and
`latents_bootstrap_failed` when fewer than two resamples produced a
usable fit, whose message carries the first reason one gave. Resamples
that fail or do not converge are dropped with a
`latents_bootstrap_dropped` warning naming how many, rather than being
silently left out of the count.

## Examples

``` r
# Two separated profiles, membership driven by `x` and by the group's class.
set.seed(5)
school <- rep(seq_len(16), each = 8)
high_class <- rep(rep(c(FALSE, TRUE), length.out = 16), each = 8)
x <- rnorm(128)
profile <- ifelse(
  runif(128) < plogis(-1 + 2 * high_class + 0.8 * x), 2L, 1L
)
example_data <- data.frame(
  school = school, x = x,
  y1 = rnorm(128, ifelse(profile == 2L, 2, -2), 0.7),
  y2 = rnorm(128, ifelse(profile == 2L, 1.5, -1.5), 0.7)
)
fit <- multilpa(example_data, c("y1", "y2"), "school", n_profiles = 2,
                n_group_classes = 2, profile_covariates = "x",
                n_starts = 2, seed = 1)
parameter_inference(fit, example_data)
#>          level       outcome          term   parameter   estimate
#> 1  measurement     profile_1            y1        mean  1.9190929
#> 2  measurement     profile_1            y2        mean  1.4517555
#> 3  measurement     profile_2            y1        mean -1.8353900
#> 4  measurement     profile_2            y2        mean -1.4838862
#> 5  measurement     profile_1            y1    variance  0.3992639
#> 6  measurement     profile_1            y2    variance  0.5601752
#> 7  measurement     profile_2            y1    variance  0.4605175
#> 8  measurement     profile_2            y2    variance  0.5202968
#> 9      profile     profile_1 group_class_1 coefficient -0.7082134
#> 10     profile     profile_1 group_class_2 coefficient  1.2934651
#> 11     profile     profile_1             x coefficient  0.7401938
#> 12       group group_class_1   (Intercept) coefficient  0.2318505
#>    standard_error   statistic       p_value    p_adjusted   conf_low
#> 1      0.07606996  25.2279984 1.975161e-140 1.975161e-140  1.7699986
#> 2      0.09010398  16.1120018  2.101105e-58  2.101105e-58  1.2751549
#> 3      0.08835082 -20.7738883  7.457039e-96  7.457039e-96 -2.0085544
#> 4      0.09390811 -15.8014691  3.039500e-56  3.039500e-56 -1.6679427
#> 5      0.06798158          NA            NA            NA  0.2660224
#> 6      0.09537546          NA            NA            NA  0.3732428
#> 7      0.08480659          NA            NA            NA  0.2942996
#> 8      0.09579579          NA            NA            NA  0.3325405
#> 9      0.33031825  -2.1440334  3.203021e-02  3.203021e-02 -1.3556252
#> 10     0.41632645   3.1068530  1.890904e-03  1.890904e-03  0.4774802
#> 11     0.22710388   3.2592740  1.116977e-03  1.116977e-03  0.2950784
#> 12     0.59791884   0.3877624  6.981918e-01  6.981918e-01 -0.9400489
#>      conf_high
#> 1   2.06818734
#> 2   1.62835603
#> 3  -1.66222558
#> 4  -1.29982964
#> 5   0.53250532
#> 6   0.74710769
#> 7   0.62673536
#> 8   0.70805306
#> 9  -0.06080148
#> 10  2.10944992
#> 11  1.18530922
#> 12  1.40374986
set.seed(42)
example_data <- data.frame(
  school = rep(seq_len(10), each = 10),
  score_a = rnorm(100), score_b = rnorm(100)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
parameter_inference(fit)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16590113
#> 2  measurement     profile_1       score_b        mean -0.07256241
#> 3  measurement     profile_2       score_a        mean -2.34985548
#> 4  measurement     profile_2       score_b        mean -0.35398815
#> 5  measurement     profile_1       score_a    variance  0.78554761
#> 6  measurement     profile_1       score_b    variance  0.71238675
#> 7  measurement     profile_2       score_a    variance  0.22457181
#> 8  measurement     profile_2       score_b    variance  2.46626299
#> 9      profile     profile_1 group_class_1 probability  0.94697964
#> 10     profile     profile_2 group_class_1 probability  0.05302036
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10676915  1.5538302 1.202249e-01 1.202249e-01 -0.04336256  0.3751648
#> 2      0.09176907 -0.7907066 4.291152e-01 4.291152e-01 -0.25242649  0.1073017
#> 3      0.38502602 -6.1031083 1.040253e-09 1.040253e-09 -3.10449262 -1.5952184
#> 4      0.79055991 -0.4477689 6.543200e-01 6.543200e-01 -1.90345711  1.1954808
#> 5      0.14679694         NA           NA           NA  0.49783091  1.0732643
#> 6      0.10792281         NA           NA           NA  0.50086192  0.9239116
#> 7      0.22279732         NA           NA           NA -0.21210290  0.6612465
#> 8      1.67675021         NA           NA           NA -0.82010704  5.7526330
#> 9      0.03531788         NA           NA           NA  0.87775787  1.0162014
#> 10     0.03531788         NA           NA           NA -0.01620141  0.1222421
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# Many tests in one table: name the correction, do not apply one by stealth.
parameter_inference(fit, adjust = "BH")
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16590113
#> 2  measurement     profile_1       score_b        mean -0.07256241
#> 3  measurement     profile_2       score_a        mean -2.34985548
#> 4  measurement     profile_2       score_b        mean -0.35398815
#> 5  measurement     profile_1       score_a    variance  0.78554761
#> 6  measurement     profile_1       score_b    variance  0.71238675
#> 7  measurement     profile_2       score_a    variance  0.22457181
#> 8  measurement     profile_2       score_b    variance  2.46626299
#> 9      profile     profile_1 group_class_1 probability  0.94697964
#> 10     profile     profile_2 group_class_1 probability  0.05302036
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10676915  1.5538302 1.202249e-01 2.404498e-01 -0.04336256  0.3751648
#> 2      0.09176907 -0.7907066 4.291152e-01 5.721536e-01 -0.25242649  0.1073017
#> 3      0.38502602 -6.1031083 1.040253e-09 4.161013e-09 -3.10449262 -1.5952184
#> 4      0.79055991 -0.4477689 6.543200e-01 6.543200e-01 -1.90345711  1.1954808
#> 5      0.14679694         NA           NA           NA  0.49783091  1.0732643
#> 6      0.10792281         NA           NA           NA  0.50086192  0.9239116
#> 7      0.22279732         NA           NA           NA -0.21210290  0.6612465
#> 8      1.67675021         NA           NA           NA -0.82010704  5.7526330
#> 9      0.03531788         NA           NA           NA  0.87775787  1.0162014
#> 10     0.03531788         NA           NA           NA -0.01620141  0.1222421
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# A structure that constrains the shape has no Wald chart. The bootstrap
# resamples schools and refits inside the same family, so it reports one.
shaped <- multilpa(example_data, c("score_a", "score_b"), "school",
                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
                   volume = "varying", shape = "equal", orientation = "axis")
parameter_inference(shaped, method = "bootstrap", iter = 10, n_starts = 1,
                    seed = 1)
#> Warning: A variance reached min_variance; this is a bound-active constrained fit.
#> Warning: A profile or group class has effective membership below one.
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.25072096
#> 2  measurement     profile_1       score_b        mean -0.06201221
#> 3  measurement     profile_2       score_a        mean -1.35056544
#> 4  measurement     profile_2       score_b        mean -0.24893253
#> 5  measurement     profile_1       score_a    variance  0.70049635
#> 6  measurement     profile_1       score_b    variance  0.70881990
#> 7  measurement     profile_2       score_a    variance  1.31170189
#> 8  measurement     profile_2       score_b    variance  1.32728800
#> 9      profile     profile_1 group_class_1 probability  0.86373072
#> 10     profile     profile_2 group_class_1 probability  0.13626928
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error statistic p_value p_adjusted     conf_low   conf_high
#> 1       0.2901831        NA      NA         NA  0.126069083  0.95194025
#> 2       0.4654784        NA      NA         NA -1.235677029  0.08222165
#> 3       0.8585431        NA      NA         NA -2.611810624 -0.11067044
#> 4       1.2097879        NA      NA         NA -1.606513935  2.16366280
#> 5       0.2542890        NA      NA         NA  0.079133222  0.78339327
#> 6       0.2229643        NA      NA         NA  0.133035001  0.78017714
#> 7       0.6034291        NA      NA         NA  0.008600376  1.73200855
#> 8       0.5990528        NA      NA         NA  0.008641422  1.68372007
#> 9       0.2976784        NA      NA         NA  0.143494035  0.98374406
#> 10      0.2976784        NA      NA         NA  0.016255939  0.85650597
#> 11      0.0000000        NA      NA         NA  1.000000000  1.00000000
```
