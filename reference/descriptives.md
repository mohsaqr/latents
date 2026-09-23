# Describe the variables a model uses, or will use

One row per variable, with the summaries you would look at before
fitting and the one that decides whether a two-level model has anything
to find: the intraclass correlation, the share of each variable's
variance that lies *between* groups rather than within them.

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

  A data frame, or a fitted model of this package. A fit carries the
  columns it was built from, so `vars` and `id` are read from it.

- ...:

  Passed between methods.

- vars:

  Character vector of column names to describe. Repeated names are
  described once. For a data frame, defaults to every numeric column
  that is not `id`.

- id:

  Optional. The column identifying the groups rows are nested in.
  Supplying it adds `n_groups` and `icc`; omitting it leaves them out.

- by:

  Optional. `"profile"` or `"group_class"` on a fitted model, or a
  column name on a data frame: describe each variable once per level of
  it, rather than once overall. Rows whose stratifier is `NA` are
  neither dropped nor folded into another level: they form their own
  stratum, labelled `NA` in the `by` column and reported last, so that
  the strata account for every input row exactly once. The stratifier is
  never compared with `==`, so a missing value cannot select rows it
  does not belong to.

## Value

A base `data.frame`, one row per variable, or one row per variable and
`by` level, with the columns

- `variable`:

  character: the column described.

- `n`:

  integer: observations with a value.

- `n_missing`:

  integer: observations without one.

- `mean`,`sd`,`min`,`max`:

  numeric, and `NA_real_` for a variable that is not numeric.

- `n_distinct`:

  integer: distinct observed values, which is the useful summary where a
  mean is not.

- `n_groups`,`icc`:

  present only when `id` is known. `icc` is the one-way random-effects
  estimate, and may be slightly negative when the between-group mean
  square falls below the within-group one; that is the estimator
  reporting no group structure, not an error. Rows with a missing group
  ID do not contribute to either statistic.

A `by` column comes first, named after what it splits on and holding
that level's own value, `NA` for the stratum of rows whose stratifier is
missing. For every variable, `n + n_missing` summed over the strata
equals the number of input rows; that invariant is asserted before the
table is returned. Describing by a class the model assigned keeps that
assignment out of the described frame, so an indicator of the caller's
own that is named `profile` or `group_class` is summarized as itself,
not replaced by the assignment. A `by` whose name is one of the
summary's own columns cannot be told apart from them and raises
`latents_bad_data`.

## Details

An ICC near zero suggests little between-group mean variation in that
variable. It does not rule out group classes that differ in other
variables, profile prevalence, or other features of their distributions.

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
