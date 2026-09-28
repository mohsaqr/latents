# Coerce a fitted multilevel latent profile model to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. The other tables are named rather than positional, so they belong
to
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
which takes `what` and refuses a name this object has not.

## Usage

``` r
# S3 method for class 'multilpa'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A fitted multilevel latent profile model.

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

A base `data.frame`: the Gaussian measurement model, one row per profile
and continuous indicator.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
for every other table this object holds.

## Examples

``` r
fit <- multilpa(
  course_engagement,
  vars = c("browse", "lectures", "forum_read", "forum_post", "attendance"),
  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
  seed = 1
)
as.data.frame(fit)
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5395774 0.5898715          0.7680309
#> 2        1   lectures  0.4274032 0.8446922          0.9190714
#> 3        1 forum_read  0.6198249 0.4814742          0.6938834
#> 4        1 forum_post  0.5304889 0.6887894          0.8299334
#> 5        1 attendance  0.6490499 0.3841262          0.6197792
#> 6        2     browse -0.7655996 0.5284686          0.7269584
#> 7        2   lectures -0.6064852 0.5383733          0.7337393
#> 8        2 forum_read -0.8792135 0.3631131          0.6025887
#> 9        2 forum_post -0.7526882 0.4211624          0.6489703
#> 10       2 attendance -0.9206598 0.3736294          0.6112523
#>    mean_standard_error variance_standard_error
#> 1           0.02710160              0.02948401
#> 2           0.03248649              0.04210053
#> 3           0.02473030              0.02439626
#> 4           0.02953027              0.03481644
#> 5           0.02226482              0.01955857
#> 6           0.03121145              0.03221701
#> 7           0.03093059              0.03209025
#> 8           0.02621510              0.02279296
#> 9           0.02768122              0.02607601
#> 10          0.02642702              0.02278996
```
