# Coerce a latent transition model summary to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. A summary carries every table the object it describes can
produce, and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
names them.

## Usage

``` r
# S3 method for class 'summary_multilpa_transitions'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `summary_multilpa_transitions`.

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
moves <- lta(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, time = "sequence",
  n_starts = 2, seed = 1
)
as.data.frame(summary(moves))
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5417725 0.5922676          0.7695892
#> 2        1   lectures  0.4284615 0.8380434          0.9154471
#> 3        1 forum_read  0.6203240 0.4829807          0.6949681
#> 4        1 forum_post  0.5313571 0.6862906          0.8284266
#> 5        1 attendance  0.6460329 0.3894969          0.6240968
#> 6        2     browse -0.7686430 0.5170201          0.7190411
#> 7        2   lectures -0.6079305 0.5447719          0.7380867
#> 8        2 forum_read -0.8798404 0.3589808          0.5991500
#> 9        2 forum_post -0.7538504 0.4216477          0.6493441
#> 10       2 attendance -0.9162947 0.3795494          0.6160758
```
