# Convert latents results to a data frame

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the primary table of a latents result: the measurement parameters of a
fitted model, the candidates of an enumeration, the test of a bootstrap,
and so on. It is the same table as `get_results(x)`;
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
with `what =` returns the others.

## Usage

``` r
# S3 method for class 'multilpa'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'multilpa_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'multilpa_start'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'multilpa_covariates'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary_multilpa_additive'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary_multilpa_covariates'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary_multilpa_enumeration'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'multilpa_bootstrap_lrt'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary_multilpa_bootstrap_lrt'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'latents_table'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary_multilpa'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'multilpa_transitions'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'summary_multilpa_transitions'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object returned by a latents function, or the
  [`summary()`](https://rdrr.io/r/base/summary.html) of one.

- row.names:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html);
  `NULL` gives default row names.

- optional:

  Ignored; present for compatibility with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- ...:

  Must be empty. An argument here raises `latents_bad_argument` naming
  it, rather than being dropped.

## Value

A base `data.frame`, one row per parameter, candidate or test depending
on the object.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
[latents-print](https://pak.dynasite.org/latents/reference/latents-print.md),
[latents-summary](https://pak.dynasite.org/latents/reference/latents-summary.md).

## Examples

``` r
fit <- multilpa(course_engagement, c("browse", "lectures", "forum_read"),
                "student", n_profiles = 2, n_group_classes = 2,
                n_starts = 2, seed = 1)
head(as.data.frame(fit))
#>   profile  indicator       mean  variance standard_deviation
#> 1       1     browse  0.5388561 0.5896074          0.7678590
#> 2       1   lectures  0.4194852 0.8397745          0.9163921
#> 3       1 forum_read  0.6130248 0.4870705          0.6979044
#> 4       2     browse -0.7837817 0.5052548          0.7108127
#> 5       2   lectures -0.6102031 0.5489694          0.7409247
#> 6       2 forum_read -0.8914110 0.3495553          0.5912320
#>   mean_standard_error variance_standard_error
#> 1          0.02838311              0.03074947
#> 2          0.03346109              0.04258343
#> 3          0.02749619              0.02753482
#> 4          0.03381727              0.03432265
#> 5          0.03327362              0.03480293
#> 6          0.02779554              0.02342192
```
