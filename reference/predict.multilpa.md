# Classify new observations with a fitted model

Evaluates the fitted mixture on new rows: the E-step of the fit, at its
estimates, applied to data it has not seen. Rows go through the fit's
own preparation – the same categorical levels, the same centering, the
same missing-data handling – so a prediction on the training data
reproduces the fit's own posteriors exactly.

## Usage

``` r
# S3 method for class 'multilpa'
predict(object, newdata = NULL, type = c("class", "posterior", "density"), ...)
```

## Arguments

- object:

  A fit from
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md)
  without membership covariates.

- newdata:

  `NULL` for the fitted rows, or a data frame with the fit's indicator
  columns (and its `id` column for a two-level fit).

- type:

  `"class"`: the modal profile and its probability, and the modal group
  class. `"posterior"`: the same with every profile's probability.
  `"density"`: each row's log density under the fitted mixture, marginal
  over the group classes (the density of a row from a new group), for
  scoring, outlier detection or held-out likelihood.

- ...:

  Unused.

## Value

A base `data.frame` with one row per row of `newdata`: `row` (its
position in `newdata`), the `id` column for a two-level fit, and

- for `"class"`: `profile` (`0` is the noise component of a
  `noise = TRUE` fit) and its `posterior`, and for two-level fits
  `group_class` and `group_posterior`;

- for `"posterior"`: those columns and one `probability_profile_k`
  column per profile (and `probability_noise`);

- for `"density"`: `log_density`.

## Groups

A two-level model classifies a row using its group: the group class is
inferred from all of that group's rows, and it shifts the row's profile
probabilities. The groups in `newdata` (the fit's `id` column) are
treated as *new* groups, each classified from its own rows alone, which
is the case of applying a model to a new cohort. A single-level fit, or
a fit given `id = NULL`, treats every row as its own unit.

## Conditions

`latents_bad_data` when `newdata` lacks an indicator or the `id` column,
has a missing value the fit cannot handle (`missing = "error"`), a
category the fit never saw, or a non-numeric continuous indicator.

## Examples

``` r
activity <- c("browse", "lectures", "forum_read")
training <- subset(course_engagement, student <= 60)
fit <- multilpa(training, activity, id = "student", n_profiles = 2,
                n_group_classes = 2, n_starts = 2, seed = 1)
new_students <- subset(course_engagement, student > 100)
head(predict(fit, new_students))
#>   row student profile posterior group_class group_posterior
#> 1   1     101       1 0.9999700           2       0.9999999
#> 2   2     101       1 0.9386987           2       0.9999999
#> 3   3     101       1 0.9855908           2       0.9999999
#> 4   4     101       1 0.9999469           2       0.9999999
#> 5   5     101       1 0.9998584           2       0.9999999
#> 6   6     101       1 0.9982339           2       0.9999999
head(predict(fit, new_students, type = "density"))
#>   row log_density
#> 1   1   -4.079463
#> 2   2   -3.333001
#> 3   3   -4.768241
#> 4   4   -3.541111
#> 5   5   -2.854424
#> 6   6   -3.886470
```
