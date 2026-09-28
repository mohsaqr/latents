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
#> 1        1     browse  0.5923990 0.5158958          0.7182588
#> 2        1   lectures  0.4581926 0.8189827          0.9049767
#> 3        1 forum_read  0.6002516 0.4758879          0.6898463
#> 4        1 forum_post  0.6050046 0.6336299          0.7960087
#> 5        1 attendance  0.6665145 0.4406607          0.6638228
#> 6        2     browse -0.7884517 0.5794487          0.7612153
#> 7        2   lectures -0.5808858 0.5564744          0.7459721
#> 8        2 forum_read -0.8371933 0.3376882          0.5811095
#> 9        2 forum_post -0.8129092 0.3733858          0.6110530
#> 10       2 attendance -0.9009123 0.3996481          0.6321772
#>    mean_standard_error variance_standard_error
#> 1           0.04431258              0.04525489
#> 2           0.05688120              0.07212904
#> 3           0.04293518              0.04169415
#> 4           0.05064070              0.05855360
#> 5           0.04268709              0.04016219
#> 6           0.07019072              0.07725965
#> 7           0.06369550              0.06634468
#> 8           0.05353799              0.04490641
#> 9           0.05381386              0.04880719
#> 10          0.05505046              0.04825216
```
