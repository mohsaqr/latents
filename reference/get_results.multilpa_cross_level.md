# Tables of a cross-level fit

Tidy tables of a fit from `multilpa(family = "restricted_cross_level")`
or `multilpa(family = "full_cross_level")`.

## Usage

``` r
# S3 method for class 'multilpa_cross_level'
get_results(x, what = "profiles", ...)

# S3 method for class 'multilpa_cross_level'
as.data.frame(x, row.names = NULL, optional = FALSE, what = "profiles", ...)

# S3 method for class 'multilpa_cross_level'
print(x, ...)

# S3 method for class 'multilpa_cross_level'
summary(object, ...)

# S3 method for class 'multilpa_cross_level'
logLik(object, ...)

# S3 method for class 'multilpa_cross_level'
nobs(object, ...)
```

## Arguments

- x:

  A `multilpa_cross_level` fit.

- what:

  Which table; `"all"` returns a named list of every table.

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

- object:

  A `multilpa_cross_level` fit.

## Value

A base `data.frame`, or a named list of them for `what = "all"`.

## Tables

- `profiles`:

  One row per individual profile and indicator: `level = "individual"`,
  `profile`, `indicator`, `mean`, `variance`; then one row per group
  class and indicator for the group means: `level = "group"`, the class
  in `profile`, and the class's `mean` and `variance` of the group
  means.

- `group_classes`:

  One row per group class: `weight`, `count` (summed posterior) and
  `n_assigned`.

- `composition`:

  One row per group class and profile: the share of the class's members
  in each profile. Model-implied for the full family; for the restricted
  family a posterior-weighted description, since its prevalences do not
  depend on the class.

- `groups`:

  One row per group: size, modal group class, its posterior and every
  class probability.

- `assignments`:

  One row per data row: its modal `profile` and posterior, and its
  group's modal `group_class` (inherited).

- `fit`:

  One row: family, log likelihood, `n_parameters`, `aic`, `bic` (by
  groups), `bic_individual` (by rows), sizes, convergence and
  `boundary`. The likelihood is the working likelihood of the manifest
  specification (the group means are computed from the same ratings), so
  it is comparable only between cross-level fits of the same data.

- `starts`:

  The start records (full family) or the two parts' likelihoods
  (restricted family).

## Examples

``` r
set.seed(1)
ratings <- data.frame(
  team = rep(seq_len(40), each = 6),
  climate = rep(rnorm(40, rep(c(-1, 1), each = 20), 0.5), each = 6) +
    rnorm(240, sample(c(-1, 1), 240, replace = TRUE), 0.6))
fit <- multilpa(ratings, "climate", "team", n_profiles = 2,
                n_group_classes = 2, family = "full_cross_level",
                n_starts = 3, seed = 1)
get_results(fit)
#>        level       profile indicator       mean  variance
#> 1 individual     profile_1   climate -1.0905752 1.3563599
#> 2 individual     profile_2   climate  1.0014259 1.1907874
#> 3      group group_class_1   climate -1.0635988 0.1439940
#> 4      group group_class_2   climate  0.9744799 0.3113962
get_results(fit, "composition")
#>     group_class   profile      share
#> 1 group_class_1 profile_1 0.98618573
#> 2 group_class_2 profile_1 0.01363281
#> 3 group_class_1 profile_2 0.01381427
#> 4 group_class_2 profile_2 0.98636719
```
