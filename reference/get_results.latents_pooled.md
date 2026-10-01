# Tables of a pooled multiply imputed fit

Tables of a pooled multiply imputed fit

## Usage

``` r
# S3 method for class 'latents_pooled'
get_results(x, what = c("estimates", "imputations", "fits"), ...)

# S3 method for class 'latents_pooled'
as.data.frame(x, row.names = NULL, optional = FALSE, what = "estimates", ...)
```

## Arguments

- x:

  A `latents_pooled` object from
  [`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md).

- what:

  `"estimates"` (the pooled table), `"imputations"` (every imputation's
  aligned estimates) or `"fits"` (one row per imputation).

- ...:

  Unused.

- row.names, optional:

  Unused; part of the generic.

## Value

A base `data.frame`; see
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
for the columns.

## Examples

``` r
set.seed(1)
# Two stand-in completed data sets. With real missing covariates, pass the
# `mids` object mice::mice() returns instead.
completed <- lapply(1:2, function(i) {
  within(subset(course_engagement, student <= 40),
         previous_grade <- previous_grade + rnorm(length(previous_grade), sd = 0.1))
})
pooled <- pool_imputations(completed, c("browse", "lectures", "forum_read"),
                           "student", n_profiles = 2, n_group_classes = 1,
                           profile_covariates = "previous_grade",
                           n_starts = 2, seed = 1)
get_results(pooled, "imputations")
#>    imputation       level   outcome           term   parameter   estimate
#> 1           1 measurement profile_1         browse        mean  0.5666028
#> 2           1 measurement profile_1       lectures        mean  0.4572288
#> 3           1 measurement profile_1     forum_read        mean  0.6304221
#> 4           1 measurement profile_2         browse        mean -0.8409613
#> 5           1 measurement profile_2       lectures        mean -0.6253489
#> 6           1 measurement profile_2     forum_read        mean -0.9166660
#> 7           1 measurement profile_1         browse    variance  0.5415729
#> 8           1 measurement profile_1       lectures    variance  0.8140909
#> 9           1 measurement profile_1     forum_read    variance  0.4746033
#> 10          1 measurement profile_2         browse    variance  0.4886298
#> 11          1 measurement profile_2       lectures    variance  0.5646644
#> 12          1 measurement profile_2     forum_read    variance  0.3389937
#> 13          1     profile profile_1    (Intercept) coefficient  0.6153411
#> 14          1     profile profile_1 previous_grade coefficient  0.6539366
#> 15          2 measurement profile_1         browse        mean  0.5680516
#> 16          2 measurement profile_1       lectures        mean  0.4583964
#> 17          2 measurement profile_1     forum_read        mean  0.6324261
#> 18          2 measurement profile_2         browse        mean -0.8384320
#> 19          2 measurement profile_2       lectures        mean -0.6235001
#> 20          2 measurement profile_2     forum_read        mean -0.9146311
#> 21          2 measurement profile_1         browse    variance  0.5411993
#> 22          2 measurement profile_1       lectures    variance  0.8134231
#> 23          2 measurement profile_1     forum_read    variance  0.4732872
#> 24          2 measurement profile_2         browse    variance  0.4893585
#> 25          2 measurement profile_2       lectures    variance  0.5664937
#> 26          2 measurement profile_2     forum_read    variance  0.3393972
#> 27          2     profile profile_1    (Intercept) coefficient  0.6130559
#> 28          2     profile profile_1 previous_grade coefficient  0.6493271
#>    standard_error
#> 1      0.04200472
#> 2      0.05093705
#> 3      0.04056609
#> 4      0.05750577
#> 5      0.05838584
#> 6      0.04884813
#> 7      0.04380687
#> 8      0.06409682
#> 9      0.03938824
#> 10     0.05737633
#> 11     0.06250250
#> 12     0.04182692
#> 13     0.10903362
#> 14     0.11061869
#> 15     0.04199082
#> 16     0.05093511
#> 17     0.04048777
#> 18     0.05732361
#> 19     0.05832402
#> 20     0.04871261
#> 21     0.04377140
#> 22     0.06409548
#> 23     0.03927968
#> 24     0.05728737
#> 25     0.06256244
#> 26     0.04171046
#> 27     0.10873868
#> 28     0.10914437
```
