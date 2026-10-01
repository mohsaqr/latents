# Tables of a latent transition fit with occasion- or covariate-dependent transitions or occasion-specific measurement

The result of
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
`transitions = "occasion"`, `transition_covariates`,
`initial_covariates` or `measurement = "occasion"`.

## Usage

``` r
# S3 method for class 'multilpa_lta'
get_results(
  x,
  what = "transitions",
  level = 0.95,
  vcov_type = c("observed", "robust", "opg"),
  ...
)

# S3 method for class 'multilpa_lta'
as.data.frame(x, row.names = NULL, optional = FALSE, what = "transitions", ...)

# S3 method for class 'multilpa_lta'
print(x, ...)

# S3 method for class 'multilpa_lta'
summary(object, ...)

# S3 method for class 'multilpa_lta'
coef(object, ...)

# S3 method for class 'multilpa_lta'
vcov(object, type = c("observed", "robust", "opg"), ...)

# S3 method for class 'multilpa_lta'
logLik(object, ...)

# S3 method for class 'multilpa_lta'
nobs(object, ...)
```

## Arguments

- x:

  A `multilpa_lta` fit.

- what:

  Which table; `"all"` for a named list of every table.

- level:

  Confidence level of the coefficient intervals.

- vcov_type:

  `"observed"`, `"robust"` (CR0 over groups) or `"opg"`.

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

- object:

  A `multilpa_lta` fit.

- type:

  Covariance type, as `vcov_type`.

## Value

A base `data.frame`, or a named list of them for `what = "all"`.

## Tables

- `transitions`:

  Model-implied transition probabilities, one row per group class,
  origin `from` and destination `to`; with occasion-varying or
  covariate-dependent transitions, one set per destination `occasion`,
  averaged over the groups observed there.

- `transition_coefficients`:

  One row per group class, move `from` -\> `to` and design `term`: the
  log odds of that move rather than staying in `from`, its standard
  error, Wald statistic, p-value, interval and `odds_ratio`. Terms are
  `(Intercept)` or `occasion_t` and each covariate column.

- `second_order_transitions`, `second_order_coefficients`:

  With `order = 2`: the probability of moving `to` given the profiles at
  the two previous occasions (`previous`, `from`), averaged over the
  groups observed from the third occasion on, and its logit coefficients
  (reference: staying in `from`). The `transitions` tables then describe
  the first move only.

- `initial`, `initial_coefficients`:

  The initial distribution (averaged over groups) and its logit
  coefficients, reference: the last profile.

- `profiles`:

  Means and variances per profile and indicator, per `occasion` when
  measurement is occasion-specific.

- `group_classes`, `assignments`, `fit`, `starts`:

  As for other fits.

## Examples

``` r
# Does the grade in the previous course shift the chance of switching
# engagement profile? `previous_grade` changes from course to course.
fit <- lta(course_engagement, c("browse", "lectures"), "student",
           n_profiles = 2, time = "sequence",
           transition_covariates = "previous_grade", n_starts = 2, seed = 1)
get_results(fit, "transition_coefficients")
#>     group_class      from        to           term    estimate standard_error
#> 1 group_class_1 profile_1 profile_2    (Intercept) -1.90765931      0.1808421
#> 2 group_class_1 profile_1 profile_2 previous_grade  0.05657504      0.1968925
#> 3 group_class_1 profile_2 profile_1    (Intercept) -1.42193554      0.1752735
#> 4 group_class_1 profile_2 profile_1 previous_grade  0.82345856      0.1996412
#>     statistic      p_value   conf_low  conf_high odds_ratio
#> 1 -10.5487573 5.147338e-26 -2.2621033 -1.5532153  0.1484274
#> 2   0.2873397 7.738522e-01 -0.3293272  0.4424773  1.0582060
#> 3  -8.1126672 4.952055e-16 -1.7654653 -1.0784058  0.2412466
#> 4   4.1246921 3.712311e-05  0.4321690  1.2147482  2.2783661
get_results(fit, "transitions")
#>      group_class occasion      from        to probability
#> 1  group_class_1        2 profile_1 profile_1   0.8704973
#> 2  group_class_1        2 profile_1 profile_2   0.1295027
#> 3  group_class_1        2 profile_2 profile_1   0.2235487
#> 4  group_class_1        2 profile_2 profile_2   0.7764513
#> 5  group_class_1        3 profile_1 profile_1   0.8701445
#> 6  group_class_1        3 profile_1 profile_2   0.1298555
#> 7  group_class_1        3 profile_2 profile_1   0.2349405
#> 8  group_class_1        3 profile_2 profile_2   0.7650595
#> 9  group_class_1        4 profile_1 profile_1   0.8710268
#> 10 group_class_1        4 profile_1 profile_2   0.1289732
#> 11 group_class_1        4 profile_2 profile_1   0.2194414
#> 12 group_class_1        4 profile_2 profile_2   0.7805586
#> 13 group_class_1        5 profile_1 profile_1   0.8704207
#> 14 group_class_1        5 profile_1 profile_2   0.1295793
#> 15 group_class_1        5 profile_2 profile_1   0.2264764
#> 16 group_class_1        5 profile_2 profile_2   0.7735236
#> 17 group_class_1        6 profile_1 profile_1   0.8706082
#> 18 group_class_1        6 profile_1 profile_2   0.1293918
#> 19 group_class_1        6 profile_2 profile_1   0.2228429
#> 20 group_class_1        6 profile_2 profile_2   0.7771571
#> 21 group_class_1        7 profile_1 profile_1   0.8690997
#> 22 group_class_1        7 profile_1 profile_2   0.1309003
#> 23 group_class_1        7 profile_2 profile_1   0.2542472
#> 24 group_class_1        7 profile_2 profile_2   0.7457528
#> 25 group_class_1        8 profile_1 profile_1   0.8716949
#> 26 group_class_1        8 profile_1 profile_2   0.1283051
#> 27 group_class_1        8 profile_2 profile_1   0.1970487
#> 28 group_class_1        8 profile_2 profile_2   0.8029513
#> 29 group_class_1        9 profile_1 profile_1   0.8705307
#> 30 group_class_1        9 profile_1 profile_2   0.1294693
#> 31 group_class_1        9 profile_2 profile_1   0.2286279
#> 32 group_class_1        9 profile_2 profile_2   0.7713721
#> 33 group_class_1       10 profile_1 profile_1   0.8709648
#> 34 group_class_1       10 profile_1 profile_2   0.1290352
#> 35 group_class_1       10 profile_2 profile_1   0.2149162
#> 36 group_class_1       10 profile_2 profile_2   0.7850838
#> 37 group_class_1       11 profile_1 profile_1   0.8702714
#> 38 group_class_1       11 profile_1 profile_2   0.1297286
#> 39 group_class_1       11 profile_2 profile_1   0.2271865
#> 40 group_class_1       11 profile_2 profile_2   0.7728135
#> 41 group_class_1       12 profile_1 profile_1   0.8700900
#> 42 group_class_1       12 profile_1 profile_2   0.1299100
#> 43 group_class_1       12 profile_2 profile_1   0.2332309
#> 44 group_class_1       12 profile_2 profile_2   0.7667691
#> 45 group_class_1       13 profile_1 profile_1   0.8704553
#> 46 group_class_1       13 profile_1 profile_2   0.1295447
#> 47 group_class_1       13 profile_2 profile_1   0.2280270
#> 48 group_class_1       13 profile_2 profile_2   0.7719730
#> 49 group_class_1       14 profile_1 profile_1   0.8718677
#> 50 group_class_1       14 profile_1 profile_2   0.1281323
#> 51 group_class_1       14 profile_2 profile_1   0.1967361
#> 52 group_class_1       14 profile_2 profile_2   0.8032639
#> 53 group_class_1       15 profile_1 profile_1   0.8707713
#> 54 group_class_1       15 profile_1 profile_2   0.1292287
#> 55 group_class_1       15 profile_2 profile_1   0.2120918
#> 56 group_class_1       15 profile_2 profile_2   0.7879082
```
