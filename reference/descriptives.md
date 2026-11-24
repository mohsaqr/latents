# Describe each variable of a data frame or a fitted model, with its variation within and between groups

`descriptives()` describes each variable of a data frame, or each
indicator of a fitted model: the number of observed and missing values,
the mean, standard deviation, minimum and maximum of a numeric variable,
and the number of distinct values. When the rows are nested in groups it
adds the number of groups, the standard deviations within and between
groups, and the intraclass correlation, the share of a variable's
variance that lies between groups.

The table can be computed for the whole data set, or once per level of
an observed grouping column (`by`). For a fitted model it can be
computed once per assigned profile or group class. Use it to inspect
variability, ranges and missingness before fitting, and to describe the
observations assigned to each class afterwards.

## Usage

``` r
descriptives(x, ...)

# S3 method for class 'data.frame'
descriptives(x, vars = NULL, id = NULL, by = NULL, ...)
```

## Arguments

- x:

  What to describe:

  - a data frame, with `vars`, `id` and `by` naming its columns;

  - a profile or class fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
    [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
    or
    [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
    with or without membership covariates;

  - an [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit
    with its default homogeneous first-order transitions, invariant
    measurement and continuous or categorical indicators.

  A fit is described on its own indicators and group identifier, so
  `vars` and `id` are not given with it. A single-level fit has no
  groups, so its table has no `n_groups`, `sd_within`, `sd_between` or
  `icc`. Anything else raises `latents_not_a_fit`, including the other
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  families (additive, dispersion and cross-level), an
  [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) fit with
  ordinal or count indicators, covariates, occasion-specific or
  second-order transitions, occasion-specific measurement, a stayer
  class or weights, and a fit from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md);
  describe the data frame it was fitted to instead.

- ...:

  Not used. For a data frame, any further argument raises a
  `latents_unused_argument` warning; for a fit, it raises
  `latents_bad_argument`.

- vars:

  For a data frame: character vector of column names to describe;
  repeated names are described once. `NULL` (default) takes every
  numeric column other than `id` and `by`. Columns must be vectors:
  matrix-valued columns raise `latents_bad_data`.

- id:

  For a data frame: optional name of the column identifying the groups
  the rows are nested in (for example the person in repeated measures).
  Supplying it adds `n_groups`, `sd_within`, `sd_between` and `icc`.

- by:

  Optional stratifier; each variable is then described once per level of
  it. `NULL` (default) describes all rows together.

  - On a data frame: the name of one of its columns.

  - On a fitted model, `"profile"`: the profile (class) each observation
    is assigned to.

  - On a fitted model, `"group_class"`: the group class the
    observation's group is assigned to; a single-level fit has one group
    class, so one stratum.

  Rows whose stratifier is `NA` form their own stratum, labelled `NA`
  and reported last, so every input row is counted in exactly one
  stratum. On a fitted model any other value raises
  `latents_bad_argument`.

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

  numeric: the mean, standard deviation, minimum and maximum of the
  observed values; `NA_real_` for a variable that is not numeric, such
  as a categorical indicator.

- `n_distinct`:

  integer: distinct observed values.

- `n_groups`:

  integer, present only when `id` is known: the groups with an observed
  value of the variable.

- `sd_within`:

  numeric, present only when `id` is known: the pooled within-group
  standard deviation, how much the values vary around their own group's
  mean (the square root of the within-group mean square). `NA` when no
  group holds more than one value, and for a variable that is not
  numeric.

- `sd_between`:

  numeric, present only when `id` is known: the standard deviation of
  the group means, each group counted once, how much the groups differ
  from each other. `NA` with fewer than two groups, and for a variable
  that is not numeric.

- `icc`:

  numeric, present only when `id` is known: the one-way random-effects
  intraclass correlation, `(MSB - MSW) / (MSB + (n0 - 1) * MSW)`, from
  the between-group and within-group mean squares of a one-way analysis
  of variance, with `n0` the average group size corrected for unequal
  sizes. It can be negative when the between-group mean square is
  smaller than the within-group mean square; this does not represent a
  negative population variance. It is `NA` for a variable that is not
  numeric, with fewer than two groups, with no group holding more than
  one value, or constant over all rows. Rows with a missing group
  identifier contribute to neither statistic.

With `by`, a first column named after the stratifier holds each
stratum's value (`NA` for the stratum of missing values). For every
variable, `n + n_missing` summed over the strata equals the number of
input rows. When describing a fit by `"profile"` or `"group_class"`, the
model's assignment is used as the stratifier without replacing an
indicator of the same name.

Errors: `latents_not_a_fit` for an `x` that is neither a data frame nor
a supported fit; `latents_bad_data` for a `vars` name the data frame
does not have, a matrix-valued column, or a `by` named like one of the
summary's own columns; `latents_nothing_to_describe` when there are no
variables to describe; `latents_bad_argument` for any other invalid
argument. Describing a fit with sampling weights warns with
`latents_unweighted_descriptives`.

## Details

Summaries by fitted class use assigned membership and do not incorporate
classification uncertainty. They describe the observations assigned to
each class and should be distinguished from the model's estimated
class-specific parameters. Rows excluded from classification are
retained in a separate stratum with a missing class label, following the
fit's `drop_empty` setting. These are unweighted descriptions of the
recorded rows, including rows with sampling weight zero. Describing a
weighted fit warns explicitly; the ICC is the unweighted one-way ANOVA
quantity, not a survey-weighted estimate.

An intraclass correlation (ICC) near zero indicates little between-group
variation in that variable's mean. It does not rule out group classes
that differ in other variables, in profile prevalence, or in other
features of the distributions.

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for model estimation and
[`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
for diagnostics of a fitted model.

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
#>   n_groups sd_within sd_between        icc
#> 1      106 0.8896664  0.5017669 0.19229259
#> 2      106 0.9448773  0.3916068 0.08779169
#> 3      106 0.8720825  0.5242171 0.22403384
#> 4      106 0.9122357  0.4582008 0.15034618
#> 5      106 0.8734499  0.5249713 0.22145308

# Missing stratification values form a separate stratum.
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
#>   n_groups sd_within sd_between        icc
#> 1      106 0.8896664  0.5017669 0.19229259
#> 2      106 0.9448773  0.3916068 0.08779169
#> 3      106 0.8720825  0.5242171 0.22403384
#> 4      106 0.9122357  0.4582008 0.15034618
#> 5      106 0.8734499  0.5249713 0.22145308
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
#>    n_groups sd_within sd_between           icc
#> 1        98 0.7583074  0.3062294  0.0162963525
#> 2        98 0.9171905  0.4167275  0.0006795179
#> 3        98 0.6925312  0.2791656 -0.0086763458
#> 4        98 0.8310652  0.3543800 -0.0073164254
#> 5        98 0.6161948  0.3061400  0.0107329239
#> 6        93 0.7255976  0.3851343 -0.0078959550
#> 7        93 0.7246577  0.4069723  0.0254729179
#> 8        93 0.5949007  0.3050038 -0.0101547583
#> 9        93 0.6531610  0.3432793 -0.0135790046
#> 10       93 0.6135168  0.3125043 -0.0088693797
```
