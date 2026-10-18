# Descriptive statistics and intraclass correlations of model variables

Returns a table with one row per variable giving the number of observed
and missing values, the mean, standard deviation, minimum, maximum and
number of distinct values, and, when groups are known, the number of
groups and the intraclass correlation (ICC): the share of the variable's
variance that lies between groups.

## Usage

``` r
descriptives(x, ...)

# S3 method for class 'data.frame'
descriptives(x, vars = NULL, id = NULL, by = NULL, ...)

# S3 method for class 'multilpa'
descriptives(x, by = NULL, ...)

# S3 method for class 'multilpa_transitions'
descriptives(x, by = NULL, ...)

# S3 method for class 'multilpa_covariates'
descriptives(x, by = NULL, ...)
```

## Arguments

- x:

  A data frame, or a fitted model of this package, whose own `vars` and
  `id` are then used. A single-level fit has no groups, so its table has
  no `n_groups` or `icc`.

- ...:

  Not used. For a data frame, any further argument raises a
  `latents_unused_argument` warning; for a fit, it raises
  `latents_bad_argument`.

- vars:

  Character vector of column names to describe; repeated names are
  described once. Columns must be vectors: matrix-valued columns raise
  `latents_bad_data`. For a data frame, defaults to every numeric column
  other than `id` and `by`.

- id:

  Optional. The column identifying the groups the rows are nested in.
  Supplying it adds `n_groups` and `icc`.

- by:

  Optional. On a fitted model, `"profile"` or `"group_class"`; on a data
  frame, a column name. Each variable is then described once per level
  of it. Rows whose stratifier is `NA` form their own stratum, labelled
  `NA` and reported last, so every input row is counted in exactly one
  stratum.

## Value

A base `data.frame`, one row per variable (or per variable and `by`
level), with the columns

- `variable`:

  character: the column described.

- `n`:

  integer: observations with a value.

- `n_missing`:

  integer: observations without one.

- `mean`,`sd`,`min`,`max`:

  numeric; `NA_real_` for a variable that is not numeric.

- `n_distinct`:

  integer: distinct observed values.

- `n_groups`,`icc`:

  present only when `id` is known. `icc` is the one-way random-effects
  estimate; it can be slightly negative when the between-group mean
  square is below the within-group one, which indicates no group
  structure. Rows with a missing group ID do not contribute to either
  statistic.

With `by`, a first column named after the stratifier holds each
stratum's value (`NA` for the stratum of missing values). For every
variable, `n + n_missing` summed over the strata equals the number of
input rows. When describing a fit by `"profile"` or `"group_class"`, the
model's assignment is used as the stratifier without replacing an
indicator of the same name. A `by` named like one of the summary's own
columns raises `latents_bad_data`. Describing no variables raises
`latents_nothing_to_describe`.

## Details

An ICC near zero indicates little between-group variation in that
variable's mean. It does not rule out group classes that differ in other
variables, in profile prevalence, or in other features of the
distributions.

## References

Bliese, P. D. (2000). Within-group agreement, non-independence, and
reliability. In K. J. Klein & S. W. J. Kozlowski (Eds.), *Multilevel
theory, research, and methods in organizations*.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
to fit,
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for the after-the-fit counterpart.

## Examples

``` r
descriptives(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student"
)
#>     variable    n n_missing          mean        sd   min  max n_distinct
#> 1     browse 1422         0 -2.812940e-05 0.9890804 -3.14 2.49        398
#> 2   lectures 1422         0 -4.219409e-05 0.9889178 -2.91 3.22        399
#> 3 forum_read 1422         0  7.032349e-05 0.9890236 -2.73 2.81        400
#> 4 forum_post 1422         0 -2.109705e-05 0.9890038 -2.61 3.10        397
#> 5 attendance 1422         0  7.735584e-05 0.9889426 -2.45 2.41        401
#>   n_groups        icc
#> 1      106 0.19229259
#> 2      106 0.08779169
#> 3      106 0.22403384
#> 4      106 0.15034618
#> 5      106 0.22145308

# A missing stratifier is a stratum of its own, not a silent loss and not
# missingness invented in a variable that has none.
descriptives(data.frame(cohort = c("A", "A", "B", NA), score = 1:4),
             vars = "score", by = "cohort")
#>   cohort variable n n_missing mean        sd min max n_distinct
#> 1      A    score 2         0  1.5 0.7071068   1   2          2
#> 2      B    score 1         0  3.0        NA   3   3          1
#> 3   <NA>    score 1         0  4.0        NA   4   4          1

fit <- multilpa(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
  seed = 1
)
descriptives(fit)
#>     variable    n n_missing          mean        sd   min  max n_distinct
#> 1     browse 1422         0 -2.812940e-05 0.9890804 -3.14 2.49        398
#> 2   lectures 1422         0 -4.219409e-05 0.9889178 -2.91 3.22        399
#> 3 forum_read 1422         0  7.032349e-05 0.9890236 -2.73 2.81        400
#> 4 forum_post 1422         0 -2.109705e-05 0.9890038 -2.61 3.10        397
#> 5 attendance 1422         0  7.735584e-05 0.9889426 -2.45 2.41        401
#>   n_groups        icc
#> 1      106 0.19229259
#> 2      106 0.08779169
#> 3      106 0.22403384
#> 4      106 0.15034618
#> 5      106 0.22145308
descriptives(fit, by = "profile")
#>    profile   variable   n n_missing       mean        sd   min  max n_distinct
#> 1        1     browse 836         0  0.5422249 0.7644929 -2.31 2.49        285
#> 2        1   lectures 836         0  0.4279785 0.9174987 -2.13 3.22        325
#> 3        1 forum_read 836         0  0.6232895 0.6895798 -1.50 2.81        273
#> 4        1 forum_post 836         0  0.5305502 0.8280756 -2.16 3.10        307
#> 5        1 attendance 836         0  0.6477632 0.6194911 -1.03 2.41        251
#> 6        2     browse 586         0 -0.7736177 0.7227891 -3.14 1.52        254
#> 7        2   lectures 586         0 -0.6106655 0.7339385 -2.91 1.43        250
#> 8        2 forum_read 586         0 -0.8890273 0.5919444 -2.73 1.15        220
#> 9        2 forum_post 586         0 -0.7569454 0.6488317 -2.61 1.59        228
#> 10       2 attendance 586         0 -0.9239249 0.6108514 -2.45 0.83        230
#>    n_groups           icc
#> 1        98  0.0162963525
#> 2        98  0.0006795179
#> 3        98 -0.0086763458
#> 4        98 -0.0073164254
#> 5        98  0.0107329239
#> 6        93 -0.0078959550
#> 7        93  0.0254729179
#> 8        93 -0.0101547583
#> 9        93 -0.0135790046
#> 10       93 -0.0088693797
```
