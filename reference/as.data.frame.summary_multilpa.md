# Coerce a multilevel LPA summary to its primary table

Plain coercion, as the base generic means it: one object, one data
frame. A summary carries every table the object it describes can
produce, and
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
names them.

## Usage

``` r
# S3 method for class 'summary_multilpa'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `summary_multilpa`.

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
  id = "student", n_profiles = 2, n_group_classes = 2, n_starts = 4,
  seed = 1
)
as.data.frame(summary(fit))
#>    profile  indicator       mean  variance standard_deviation
#> 1        1     browse  0.5395948 0.5898645          0.7680263
#> 2        1   lectures  0.4274286 0.8446803          0.9190649
#> 3        1 forum_read  0.6198475 0.4814627          0.6938751
#> 4        1 forum_post  0.5305152 0.6887698          0.8299216
#> 5        1 attendance  0.6490757 0.3841104          0.6197664
#> 6        2     browse -0.7655654 0.5284937          0.7269757
#> 7        2   lectures -0.6064745 0.5383778          0.7337424
#> 8        2 forum_read -0.8791779 0.3631402          0.6026112
#> 9        2 forum_post -0.7526676 0.4211809          0.6489846
#> 10       2 attendance -0.9206255 0.3736487          0.6112681
get_results(summary(fit), what = "counts")
#>         level class effective_count effective_proportion
#> 1 individuals     1       834.08182            0.5865554
#> 2 individuals     2       587.91818            0.4134446
#> 3      groups     1        34.47767            0.3252610
#> 4      groups     2        71.52233            0.6747390
```
