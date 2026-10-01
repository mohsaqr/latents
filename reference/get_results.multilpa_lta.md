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
#>     group_class      from        to           term   estimate standard_error
#> 1 group_class_1 profile_1 profile_2    (Intercept) -1.9079097      0.1809288
#> 2 group_class_1 profile_1 profile_2 previous_grade  0.0579463      0.1968674
#> 3 group_class_1 profile_2 profile_1    (Intercept) -1.4212912      0.1752162
#> 4 group_class_1 profile_2 profile_1 previous_grade  0.8240685      0.1995587
#>     statistic      p_value   conf_low  conf_high odds_ratio
#> 1 -10.5450838 5.352515e-26 -2.2625237 -1.5532957  0.1483902
#> 2   0.2943418 7.684968e-01 -0.3279068  0.4437994  1.0596581
#> 3  -8.1116425 4.994005e-16 -1.7647086 -1.0778737  0.2414021
#> 4   4.1294533 3.636269e-05  0.4329406  1.2151965  2.2797562
get_results(fit, "transitions")
#>      group_class occasion      from        to probability
#> 1  group_class_1        2 profile_1 profile_1   0.8705162
#> 2  group_class_1        2 profile_1 profile_2   0.1294838
#> 3  group_class_1        2 profile_2 profile_1   0.2236869
#> 4  group_class_1        2 profile_2 profile_2   0.7763131
#> 5  group_class_1        3 profile_1 profile_1   0.8701543
#> 6  group_class_1        3 profile_1 profile_2   0.1298457
#> 7  group_class_1        3 profile_2 profile_1   0.2350901
#> 8  group_class_1        3 profile_2 profile_2   0.7649099
#> 9  group_class_1        4 profile_1 profile_1   0.8710576
#> 10 group_class_1        4 profile_1 profile_2   0.1289424
#> 11 group_class_1        4 profile_2 profile_1   0.2195794
#> 12 group_class_1        4 profile_2 profile_2   0.7804206
#> 13 group_class_1        5 profile_1 profile_1   0.8704375
#> 14 group_class_1        5 profile_1 profile_2   0.1295625
#> 15 group_class_1        5 profile_2 profile_1   0.2266158
#> 16 group_class_1        5 profile_2 profile_2   0.7733842
#> 17 group_class_1        6 profile_1 profile_1   0.8706294
#> 18 group_class_1        6 profile_1 profile_2   0.1293706
#> 19 group_class_1        6 profile_2 profile_1   0.2229776
#> 20 group_class_1        6 profile_2 profile_2   0.7770224
#> 21 group_class_1        7 profile_1 profile_1   0.8690846
#> 22 group_class_1        7 profile_1 profile_2   0.1309154
#> 23 group_class_1        7 profile_2 profile_1   0.2544163
#> 24 group_class_1        7 profile_2 profile_2   0.7455837
#> 25 group_class_1        8 profile_1 profile_1   0.8717429
#> 26 group_class_1        8 profile_1 profile_2   0.1282571
#> 27 group_class_1        8 profile_2 profile_1   0.1971595
#> 28 group_class_1        8 profile_2 profile_2   0.8028405
#> 29 group_class_1        9 profile_1 profile_1   0.8705497
#> 30 group_class_1        9 profile_1 profile_2   0.1294503
#> 31 group_class_1        9 profile_2 profile_1   0.2287724
#> 32 group_class_1        9 profile_2 profile_2   0.7712276
#> 33 group_class_1       10 profile_1 profile_1   0.8709948
#> 34 group_class_1       10 profile_1 profile_2   0.1290052
#> 35 group_class_1       10 profile_2 profile_1   0.2150448
#> 36 group_class_1       10 profile_2 profile_2   0.7849552
#> 37 group_class_1       11 profile_1 profile_1   0.8702850
#> 38 group_class_1       11 profile_1 profile_2   0.1297150
#> 39 group_class_1       11 profile_2 profile_1   0.2273276
#> 40 group_class_1       11 profile_2 profile_2   0.7726724
#> 41 group_class_1       12 profile_1 profile_1   0.8700988
#> 42 group_class_1       12 profile_1 profile_2   0.1299012
#> 43 group_class_1       12 profile_2 profile_1   0.2333779
#> 44 group_class_1       12 profile_2 profile_2   0.7666221
#> 45 group_class_1       13 profile_1 profile_1   0.8704727
#> 46 group_class_1       13 profile_1 profile_2   0.1295273
#> 47 group_class_1       13 profile_2 profile_1   0.2281709
#> 48 group_class_1       13 profile_2 profile_2   0.7718291
#> 49 group_class_1       14 profile_1 profile_1   0.8719194
#> 50 group_class_1       14 profile_1 profile_2   0.1280806
#> 51 group_class_1       14 profile_2 profile_1   0.1968471
#> 52 group_class_1       14 profile_2 profile_2   0.8031529
#> 53 group_class_1       15 profile_1 profile_1   0.8707974
#> 54 group_class_1       15 profile_1 profile_2   0.1292026
#> 55 group_class_1       15 profile_2 profile_1   0.2122120
#> 56 group_class_1       15 profile_2 profile_2   0.7877880
```
