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
  subset(course_engagement, student <= 30),
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2,
  profile_covariates = "sequence", n_starts = 1, seed = 1
)
as.data.frame(summary(fit))
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse -0.7884408 0.5794605          0.7612230
#> 2        1   lectures -0.5808858 0.5564737          0.7459717
#> 3        1 forum_read -0.8371850 0.3376948          0.5811152
#> 4        1 forum_post -0.8129060 0.3733900          0.6110565
#> 5        1 attendance -0.9009078 0.3996484          0.6321775
#> 6        2     browse  0.5923999 0.5158958          0.7182589
#> 7        2   lectures  0.4581977 0.8189790          0.9049746
#> 8        2 forum_read  0.6002542 0.4758876          0.6898461
#> 9        2 forum_post  0.6050098 0.6336240          0.7960050
#> 10       2 attendance  0.6665197 0.4406563          0.6638195
```
