# Coerce a covariate-model summary to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. A summary carries every table the object it describes can
produce, and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
names them.

## Usage

``` r
# S3 method for class 'summary_multilpa_covariates'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `summary_multilpa_covariates`.

- row.names:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html);
  `NULL` gives default row names.

- optional:

  Ignored, present for generic compatibility.

- ...:

  Must be empty. An argument here raises `latents_bad_argument` naming
  it, rather than being dropped.

## Value

A base `data.frame`: one row per profile and continuous indicator.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for every other table this summary holds.

## Examples

``` r
fit <- multilpa(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2,
  profile_covariates = "sequence", n_starts = 4, seed = 1
)
as.data.frame(summary(fit))
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5395857 0.5898775          0.7680348
#> 2        1   lectures  0.4274301 0.8446567          0.9190521
#> 3        1 forum_read  0.6198229 0.4814994          0.6939016
#> 4        1 forum_post  0.5304968 0.6887881          0.8299326
#> 5        1 attendance  0.6490624 0.3841107          0.6197666
#> 6        2     browse -0.7655862 0.5284633          0.7269548
#> 7        2   lectures -0.6065034 0.5383714          0.7337380
#> 8        2 forum_read -0.8791818 0.3631316          0.6026040
#> 9        2 forum_post -0.7526746 0.4211726          0.6489781
#> 10       2 attendance -0.9206473 0.3736434          0.6112638
```
