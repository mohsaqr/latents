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
# S3 method for class 'multilpa_additive'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  ...
)

# S3 method for class 'multilpa_covariates'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

# S3 method for class 'multilpa_cross_level'
parameter_inference(x, ...)

parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

# S3 method for class 'multilpa'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  adjust = .multilpa_p_adjust_methods,
  method = c("wald", "bootstrap"),
  iter = 199L,
  n_starts = 10L,
  max_iter = 1000L,
  tol = 1e-08,
  seed = NULL,
  boundary = c("error", "fix")
)

# S3 method for class 'multilpa_lta'
parameter_inference(
  x,
  data = NULL,
  level = 0.95,
  step = 1e-04,
  vcov_type = c("observed", "robust", "opg"),
  ...
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
  groups are independent. `"opg"` inverts the outer product of the
  per-group scores (the BHHH or outer-product-of-gradients estimate),
  which estimates the information by the variance of the score rather
  than by the curvature of the likelihood. It is the estimator `glca`
  reports, so it reproduces that package's standard errors; the two
  agree with `"observed"` in large samples and can differ noticeably in
  small ones, where `"observed"` is usually the more reliable.
  `"robust"` and `"opg"` need more groups than parameters.

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
  for all fourteen covariance structures. EEI, VVI, EEE and VVV are
  estimated in log variances or log-Cholesky coordinates; the ten
  structures that constrain the volume, the shape or the orientation
  across profiles are estimated in log volumes, log shapes whose
  determinant-one constraint is written into the coordinates, and
  Cayley-transform orientations (or, for VEE and EVV, determinant-one
  log-Cholesky factors), so the information has one row per free
  parameter of the structure. The reported variances and covariances are
  then carried from those coordinates by the delta method, and their
  natural-scale covariance is singular, as it is for the class
  probabilities. `"bootstrap"` resamples *groups* with replacement,
  refits inside the same covariance family, undoes the relabelling each
  refit comes back with, and reports the percentile interval and the
  standard deviation of the replicates. It is available for all fourteen
  structures. Groups are the resampling unit rather than rows, so the
  interval carries the same independence assumption as
  `vcov_type = "robust"` and not the stronger one `"observed"` makes.

- ...:

  Unused by most methods; accepted so every method shares the generic's
  signature.

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

- boundary:

  What to do when categorical response probabilities sit on their lower
  bound, `min_probability`. There the likelihood is flat in the logit
  and the information is singular, so `"error"`, the default, refuses
  with `latents_boundary_fit`. `"fix"` holds those probabilities at the
  bound and reports every other parameter conditionally on them, the
  convention Mplus follows; the held probabilities are left out of the
  table (or, for a covariate fit, reported with `NA` standard errors)
  and named in the result's `fixed_at_bound` attribute.

## Value

A base `data.frame` with one row per free parameter and the same
columns, in the same order, that `parameter_inference()` returns for
every other fitted class: `level` (`"measurement"`, `"profile"` or
`"group"`), `outcome` (which profile or group class the coefficient
predicts, or which profile a measurement parameter belongs to), `term`,
`parameter` (`"mean"`, `"variance"`, `"covariance"`, `"response"` or
`"coefficient"`), `estimate`, `standard_error`, `statistic`, `p_value`,
`p_adjusted`, `conf_low` and `conf_high`. Variances and covariances are
reported in their natural units, with standard errors carried through
the delta method from the log and log-Cholesky coordinates on which they
are estimated. Categorical indicators contribute one `"response"` row
per profile, indicator and category, named `indicator:category`: the
probability of that category, with its standard error carried from the
multinomial logits it is estimated on, and no Wald test, as for the
covariate-free model.

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
natural covariance singular. Wald intervals for probabilities are formed
on the logit scale and for variances (and covariance diagonals) on the
log scale, then mapped back with the same standard error, so they stay
inside the parameter's range and widen towards the far side of a bound;
every other interval is `estimate +/- z * standard_error`.
[`confint()`](https://rdrr.io/r/stats/confint.html) returns the same
intervals. A probability at or below 1e-6 is treated as on its bound
(refused, or held by `boundary = "fix"`).

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
for what that conditioning means. For a
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
result, `method = "bootstrap"` does propagate it: each resample refits
both stages, and the table has a row for every parameter, the
measurement included (attribute `stages_resampled` is 2).

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
`latents_unconverged` warning rather than refused. A fit made with
`multilpa(prior = )` is refused by the Wald path with
`latents_unsupported_inference` (it sits at a posterior mode, where the
likelihood's score is not zero), and a fit with a noise component by
either path with `latents_unsupported_noise`.

`method = "bootstrap"` adds `latents_unsupported_inference` for a fit
that holds a measurement block through `multilpa(fixed = )` — the held
values came from another fit and resampling these data does not resample
them; a
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
result is bootstrapped by staging again instead — and
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
#> 9      profile     profile_1 group_class_1 coefficient -0.7082122
#> 10     profile     profile_1 group_class_2 coefficient  1.2934665
#> 11     profile     profile_1             x coefficient  0.7401937
#> 12       group group_class_1   (Intercept) coefficient  0.2318527
#>    standard_error   statistic       p_value    p_adjusted   conf_low
#> 1      0.07606996  25.2279984 1.975161e-140 1.975161e-140  1.7699986
#> 2      0.09010398  16.1120018  2.101105e-58  2.101105e-58  1.2751549
#> 3      0.08835082 -20.7738883  7.457039e-96  7.457039e-96 -2.0085544
#> 4      0.09390811 -15.8014691  3.039500e-56  3.039500e-56 -1.6679427
#> 5      0.06798158          NA            NA            NA  0.2859751
#> 6      0.09537546          NA            NA            NA  0.4012345
#> 7      0.08480659          NA            NA            NA  0.3209912
#> 8      0.09579579          NA            NA            NA  0.3626853
#> 9      0.33031825  -2.1440300  3.203049e-02  3.203049e-02 -1.3556241
#> 10     0.41632648   3.1068563  1.890883e-03  1.890883e-03  0.4774816
#> 11     0.22710387   3.2592737  1.116978e-03  1.116978e-03  0.2950783
#> 12     0.59791891   0.3877661  6.981892e-01  6.981892e-01 -0.9400469
#>      conf_high
#> 1   2.06818734
#> 2   1.62835603
#> 3  -1.66222558
#> 4  -1.29982964
#> 5   0.55743188
#> 6   0.78207700
#> 7   0.66069214
#> 8   0.74640116
#> 9  -0.06080035
#> 10  2.10945144
#> 11  1.18530911
#> 12  1.40375220
set.seed(42)
example_data <- data.frame(
  school = rep(seq_len(10), each = 10),
  score_a = rnorm(100), score_b = rnorm(100)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
parameter_inference(fit)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16569033
#> 2  measurement     profile_1       score_b        mean -0.07267264
#> 3  measurement     profile_2       score_a        mean -2.35114132
#> 4  measurement     profile_2       score_b        mean -0.35258109
#> 5  measurement     profile_1       score_a    variance  0.78587057
#> 6  measurement     profile_1       score_b    variance  0.71248993
#> 7  measurement     profile_2       score_a    variance  0.22391803
#> 8  measurement     profile_2       score_b    variance  2.46889204
#> 9      profile     profile_1 group_class_1 probability  0.94708605
#> 10     profile     profile_2 group_class_1 probability  0.05291395
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10691280  1.5497708 1.211965e-01 1.211965e-01 -0.04385490  0.3752356
#> 2      0.09184042 -0.7912926 4.287733e-01 4.287733e-01 -0.25267656  0.1073313
#> 3      0.38614364 -6.0887740 1.137786e-09 1.137786e-09 -3.10796894 -1.5943137
#> 4      0.79370695 -0.4442207 6.568830e-01 6.568830e-01 -1.90821812  1.2030559
#> 5      0.14710614         NA           NA           NA  0.54452229  1.1341915
#> 6      0.10799758         NA           NA           NA  0.52936541  0.9589631
#> 7      0.22273970         NA           NA           NA  0.03186888  1.5732992
#> 8      1.68065393         NA           NA           NA  0.65022241  9.3743738
#> 9      0.03538773         NA           NA           NA  0.81767987  0.9861939
#> 10     0.03538773         NA           NA           NA  0.01380615  0.1823201
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# Many tests in one table: name the correction, do not apply one by stealth.
parameter_inference(fit, adjust = "BH")
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.16569033
#> 2  measurement     profile_1       score_b        mean -0.07267264
#> 3  measurement     profile_2       score_a        mean -2.35114132
#> 4  measurement     profile_2       score_b        mean -0.35258109
#> 5  measurement     profile_1       score_a    variance  0.78587057
#> 6  measurement     profile_1       score_b    variance  0.71248993
#> 7  measurement     profile_2       score_a    variance  0.22391803
#> 8  measurement     profile_2       score_b    variance  2.46889204
#> 9      profile     profile_1 group_class_1 probability  0.94708605
#> 10     profile     profile_2 group_class_1 probability  0.05291395
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic      p_value   p_adjusted    conf_low  conf_high
#> 1      0.10691280  1.5497708 1.211965e-01 2.423931e-01 -0.04385490  0.3752356
#> 2      0.09184042 -0.7912926 4.287733e-01 5.716977e-01 -0.25267656  0.1073313
#> 3      0.38614364 -6.0887740 1.137786e-09 4.551145e-09 -3.10796894 -1.5943137
#> 4      0.79370695 -0.4442207 6.568830e-01 6.568830e-01 -1.90821812  1.2030559
#> 5      0.14710614         NA           NA           NA  0.54452229  1.1341915
#> 6      0.10799758         NA           NA           NA  0.52936541  0.9589631
#> 7      0.22273970         NA           NA           NA  0.03186888  1.5732992
#> 8      1.68065393         NA           NA           NA  0.65022241  9.3743738
#> 9      0.03538773         NA           NA           NA  0.81767987  0.9861939
#> 10     0.03538773         NA           NA           NA  0.01380615  0.1823201
#> 11     0.00000000         NA           NA           NA  1.00000000  1.0000000

# A structure that constrains the shape across profiles (VEI) is charted in
# its own free coordinates, so it has Wald standard errors too.
shaped <- multilpa(example_data, c("score_a", "score_b"), "school",
                   n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1,
                   volume = "varying", shape = "equal", orientation = "axis")
parameter_inference(shaped)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.25022879
#> 2  measurement     profile_1       score_b        mean -0.06203038
#> 3  measurement     profile_2       score_a        mean -1.35411736
#> 4  measurement     profile_2       score_b        mean -0.24959732
#> 5  measurement     profile_1       score_a    variance  0.70077548
#> 6  measurement     profile_1       score_b    variance  0.70904425
#> 7  measurement     profile_2       score_a    variance  1.31232456
#> 8  measurement     profile_2       score_b    variance  1.32780927
#> 9      profile     profile_1 group_class_1 probability  0.86429738
#> 10     profile     profile_2 group_class_1 probability  0.13570262
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error  statistic    p_value p_adjusted    conf_low  conf_high
#> 1      0.13357063  1.8733818 0.06101568 0.06101568 -0.01156483 0.51202242
#> 2      0.09872161 -0.6283364 0.52978357 0.52978357 -0.25552118 0.13146041
#> 3      0.72483385 -1.8681762 0.06173751 0.06173751 -2.77476561 0.06653089
#> 4      0.41275474 -0.6047110 0.54537106 0.54537106 -1.05858176 0.55938711
#> 5      0.13997195         NA         NA         NA  0.47376290 1.03656549
#> 6      0.11795100         NA         NA         NA  0.51176924 0.98236413
#> 7      0.61223539         NA         NA         NA  0.52593539 3.27453859
#> 8      0.48454076         NA         NA         NA  0.64940766 2.71490090
#> 9      0.10218086         NA         NA         NA  0.53592034 0.97231994
#> 10     0.10218086         NA         NA         NA  0.02768006 0.46407966
#> 11     0.00000000         NA         NA         NA  1.00000000 1.00000000

# The bootstrap resamples schools and refits inside the same family.
parameter_inference(shaped, method = "bootstrap", iter = 10, n_starts = 1,
                    seed = 1)
#>          level       outcome          term   parameter    estimate
#> 1  measurement     profile_1       score_a        mean  0.25022879
#> 2  measurement     profile_1       score_b        mean -0.06203038
#> 3  measurement     profile_2       score_a        mean -1.35411736
#> 4  measurement     profile_2       score_b        mean -0.24959732
#> 5  measurement     profile_1       score_a    variance  0.70077548
#> 6  measurement     profile_1       score_b    variance  0.70904425
#> 7  measurement     profile_2       score_a    variance  1.31232456
#> 8  measurement     profile_2       score_b    variance  1.32780927
#> 9      profile     profile_1 group_class_1 probability  0.86429738
#> 10     profile     profile_2 group_class_1 probability  0.13570262
#> 11       group group_class_1          <NA> probability  1.00000000
#>    standard_error statistic p_value p_adjusted    conf_low conf_high
#> 1       0.3673629        NA      NA         NA  0.03448968 1.0490512
#> 2       0.4769439        NA      NA         NA -0.30194200 1.1224914
#> 3       1.0307447        NA      NA         NA -2.68649048 0.1362221
#> 4       0.5430616        NA      NA         NA -0.69622362 0.9697939
#> 5       0.2765594        NA      NA         NA  0.20977054 1.0503136
#> 6       0.2198624        NA      NA         NA  0.15771978 0.7749002
#> 7       0.5999388        NA      NA         NA  0.03880418 1.7155546
#> 8       0.5841325        NA      NA         NA  0.04243513 1.6601809
#> 9       0.2848561        NA      NA         NA  0.16211188 0.9682255
#> 10      0.2848561        NA      NA         NA  0.03177455 0.8378881
#> 11      0.0000000        NA      NA         NA  1.00000000 1.0000000
```
