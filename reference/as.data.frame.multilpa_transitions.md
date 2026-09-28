# Coerce a fitted latent transition model to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. The other tables are named rather than positional, so they belong
to
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
which takes `what` and refuses a name this object has not.

## Usage

``` r
# S3 method for class 'multilpa_transitions'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `multilpa_transitions`.

- row.names:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html);
  `NULL` gives default row names.

- optional:

  Ignored, present for generic compatibility.

- ...:

  Must be empty. An argument here raises `latents_bad_argument` naming
  it, rather than being dropped, because `what =` used to live on this
  generic and silently returning the primary table instead of the one
  that was asked for is the one outcome worth refusing.

## Value

A base `data.frame`: one row per group class and ordered pair of
profiles.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for every other table this object holds.

## Examples

``` r
moves <- lta(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, time = "sequence",
  n_starts = 2, seed = 1
)
as.data.frame(moves)
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse -0.7686436 0.5170195          0.7190407
#> 2        1   lectures -0.6079308 0.5447709          0.7380860
#> 3        1 forum_read -0.8798410 0.3589804          0.5991497
#> 4        1 forum_post -0.7538504 0.4216477          0.6493441
#> 5        1 attendance -0.9162955 0.3795487          0.6160753
#> 6        2     browse  0.5417724 0.5922675          0.7695892
#> 7        2   lectures  0.4284612 0.8380440          0.9154474
#> 8        2 forum_read  0.6203238 0.4829806          0.6949680
#> 9        2 forum_post  0.5313566 0.6862912          0.8284269
#> 10       2 attendance  0.6460327 0.3894967          0.6240967
#>    mean_standard_error variance_standard_error
#> 1           0.03037963              0.03108758
#> 2           0.03088464              0.03254936
#> 3           0.02534797              0.02162357
#> 4           0.02739956              0.02570203
#> 5           0.02604185              0.02265529
#> 6           0.02691369              0.02930066
#> 7           0.03203979              0.04145235
#> 8           0.02452562              0.02440736
#> 9           0.02905255              0.03408914
#> 10          0.02205632              0.01961044
```
