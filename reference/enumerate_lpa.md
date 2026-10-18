# Compare latent profile models

Fits single-level latent profile models with each number of profiles in
`n_profiles` under each covariance structure in `model`, and returns
them in one table for comparison by information criteria, entropy and
the diagnostics of each fit. It is
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with `id = NULL`, and it returns the same object.

## Usage

``` r
enumerate_lpa(data, vars, n_profiles = 1:4, model = "basic", ...)
```

## Arguments

- data:

  A data frame with one row per observation.

- vars:

  Names of the indicator columns.

- n_profiles:

  Positive whole numbers of profiles to compare.

- model:

  The covariance structures: `"basic"`, the default, for the four that
  combine equal or varying variances with covariances absent or present
  (`EEI`, `VVI`, `EEE`, `VVV`); `"all"` for all 14; or their
  three-letter codes. See
  [`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md).

- ...:

  Further arguments for
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
  such as `missing`, `n_starts`, `seed`, `boundary` and `categorical`,
  except `id` and `n_group_classes`.

## Value

An object of class `multilpa_enumeration`, as
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
returns: [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
gives one row per candidate,
[`summary()`](https://rdrr.io/r/base/summary.html) the candidate each
criterion prefers,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the criteria
against the number of profiles, and
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
one fitted candidate.

## Details

For observations nested in groups, use
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with the grouping column as `id`, which also crosses the number of group
classes.

## Conditions

`latents_bad_argument` when `id` or `n_group_classes` is passed, or when
`model` names an unknown structure.

## References

Celeux, G., & Govaert, G. (1995). Gaussian parsimonious clustering
models. *Pattern Recognition*, 28(5), 781–793.

## See also

[`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) to fit one
model;
[`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md)
for categorical indicators;
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
to test one number of profiles against another.

## Examples

``` r
models <- enumerate_lpa(srl, c("cognitive_strategies", "intrinsic_value",
                               "self_efficacy", "self_regulation",
                               "test_anxiety"),
                        n_profiles = 1:3, n_starts = 2, seed = 1)
summary(models)
#> Class enumeration: 12 candidates, 12 converged, 0 failed to fit
#> 2 distinct candidate(s) are minimal under some criterion.
#> No candidate is selected automatically. Choose one convention and keep it.
#> The candidates table shows 14 of its 35 columns; get_results(x, what = "candidates") returns all of them.
#> 
#> -- candidates ------------------------------------------------------
#>      candidate n_profiles model log_likelihood n_parameters  aic  bic
#>  profiles1_EEI          1   EEI        -1247.3           10 2515 2552
#>  profiles2_EEI          2   EEI         -906.7           16 1845 1905
#>  profiles3_EEI          3   EEI         -757.7           22 1559 1641
#>  profiles1_VVI          1   VVI        -1247.3           10 2515 2552
#>  profiles2_VVI          2   VVI         -889.3           21 1821 1898
#>  profiles3_VVI          3   VVI         -735.7           32 1535 1654
#>  profiles1_EEE          1   EEE         -615.7           20 1271 1345
#>  profiles2_EEE          2   EEE         -600.9           26 1254 1350
#>  profiles3_EEE          3   EEE         -594.4           32 1253 1371
#>  profiles1_VVV          1   VVV         -615.7           20 1271 1345
#>  icl_individual profile_entropy avepp_min avepp_mean avepp_below converged
#>            2552              NA        NA         NA          NA      TRUE
#>            1949          0.8939    0.9657     0.9712           0      TRUE
#>            1707          0.9003    0.9452     0.9535           0      TRUE
#>            2552              NA        NA         NA          NA      TRUE
#>            1940          0.8991    0.9686     0.9740           0      TRUE
#>            1714          0.9089    0.9356     0.9591           0      TRUE
#>            1345              NA        NA         NA          NA      TRUE
#>            1457          0.7440    0.8801     0.9169           0      TRUE
#>            1606          0.6442    0.7797     0.8391           1      TRUE
#>            1345              NA        NA         NA          NA      TRUE
#>  boundary
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>     FALSE
#>    ... 2 more rows.  get_results(x, what = "candidates")
#> 
#> -- best ------------------------------------------------------------
#>  criterion  convention     candidate n_profiles n_group_classes model value
#>        aic        <NA> profiles3_VVV          3               1   VVV  1193
#>        kic        <NA> profiles3_VVV          3               1   VVV  1258
#>        bic      groups profiles1_EEE          1               1   EEE  1345
#>        bic individuals profiles1_EEE          1               1   EEE  1345
#>      sabic      groups profiles3_VVV          3               1   VVV  1226
#>      sabic individuals profiles3_VVV          3               1   VVV  1226
#>       caic      groups profiles1_EEE          1               1   EEE  1365
#>       caic individuals profiles1_EEE          1               1   EEE  1365
#>        awe      groups profiles1_EEE          1               1   EEE  1520
#>        awe individuals profiles1_EEE          1               1   EEE  1520
#>    ... 4 more rows.  get_results(x, what = "best")
#> 
#> -- classification --------------------------------------------------
#>      candidate n_profiles n_group_classes model       level class n_modal
#>  profiles2_EEI          2               1   EEI individuals     1     131
#>  profiles2_EEI          2               1   EEI individuals     2     169
#>  profiles3_EEI          3               1   EEI individuals     1     123
#>  profiles3_EEI          3               1   EEI individuals     2     120
#>  profiles3_EEI          3               1   EEI individuals     3      57
#>  profiles2_VVI          2               1   VVI individuals     1     140
#>  profiles2_VVI          2               1   VVI individuals     2     160
#>  profiles3_VVI          3               1   VVI individuals     1     118
#>  profiles3_VVI          3               1   VVI individuals     2     122
#>  profiles3_VVI          3               1   VVI individuals     3      60
#>  proportion_modal estimated_n estimated_proportion average_posterior
#>            0.4367      130.43               0.4348            0.9657
#>            0.5633      169.57               0.5652            0.9768
#>            0.4100      123.06               0.4102            0.9452
#>            0.4000      120.59               0.4020            0.9687
#>            0.1900       56.34               0.1878            0.9466
#>            0.4667      142.13               0.4738            0.9793
#>            0.5333      157.87               0.5262            0.9686
#>            0.3933      115.46               0.3849            0.9356
#>            0.4067      124.18               0.4139            0.9748
#>            0.2000       60.36               0.2012            0.9668
#>  odds_correct_classification relative_entropy
#>                        36.56           0.8939
#>                        32.32           0.8939
#>                        24.81           0.9003
#>                        46.01           0.9003
#>                        76.70           0.9003
#>                        52.62           0.8991
#>                        27.77           0.8991
#>                        23.22           0.9089
#>                        54.87           0.9089
#>                       115.74           0.9089
#>    ... 10 more rows.  get_results(x, what = "classification")
#> 
#> -- deltas ----------------------------------------------------------
#>      candidate n_profiles n_group_classes model  criterion value delta_profile
#>  profiles1_EEI          1               1   EEI bic_groups  2552            NA
#>  profiles2_EEI          2               1   EEI bic_groups  1905       646.836
#>  profiles3_EEI          3               1   EEI bic_groups  1641       263.844
#>  profiles1_VVI          1               1   VVI bic_groups  2552            NA
#>  profiles2_VVI          2               1   VVI bic_groups  1898       653.179
#>  profiles3_VVI          3               1   VVI bic_groups  1654       244.369
#>  profiles1_EEE          1               1   EEE bic_groups  1345            NA
#>  profiles2_EEE          2               1   EEE bic_groups  1350        -4.704
#>  profiles3_EEE          3               1   EEE bic_groups  1371       -21.098
#>  profiles1_VVV          1               1   VVV bic_groups  1345            NA
#>   from_profile delta_group_class from_group_class
#>           <NA>                NA             <NA>
#>  profiles1_EEI                NA             <NA>
#>  profiles2_EEI                NA             <NA>
#>           <NA>                NA             <NA>
#>  profiles1_VVI                NA             <NA>
#>  profiles2_VVI                NA             <NA>
#>           <NA>                NA             <NA>
#>  profiles1_EEE                NA             <NA>
#>  profiles2_EEE                NA             <NA>
#>           <NA>                NA             <NA>
#>    ... 2 more rows.  get_results(x, what = "deltas")
#> 
#> -- knees -----------------------------------------------------------
#>   direction                line  criterion          from            to n_points
#>    profiles  1 group class, EEE bic_groups profiles1_EEE profiles3_EEE        3
#>    profiles  1 group class, EEI bic_groups profiles1_EEI profiles3_EEI        3
#>    profiles  1 group class, VVI bic_groups profiles1_VVI profiles3_VVI        3
#>    profiles  1 group class, VVV bic_groups profiles1_VVV profiles3_VVV        3
#>  parameters all candidates, EEE bic_groups profiles1_EEE profiles3_EEE        3
#>  parameters all candidates, EEI bic_groups profiles1_EEI profiles3_EEI        3
#>  parameters all candidates, VVI bic_groups profiles1_VVI profiles3_VVI        3
#>  parameters all candidates, VVV bic_groups profiles1_VVV profiles3_VVV        3
#>  knee_candidate distance knee_without_largest stable
#>            <NA>       NA                 <NA>     NA
#>   profiles2_EEI   0.1487                 <NA>     NA
#>   profiles2_VVI   0.1610                 <NA>     NA
#>            <NA>       NA                 <NA>     NA
#>            <NA>       NA                 <NA>     NA
#>   profiles2_EEI   0.1487                 <NA>     NA
#>   profiles2_VVI   0.1610                 <NA>     NA
#>            <NA>       NA                 <NA>     NA
#>                                         note
#>  the criterion does not fall along this line
#>     three candidates: stability not assessed
#>     three candidates: stability not assessed
#>  the criterion does not fall along this line
#>  the criterion does not fall along this line
#>     three candidates: stability not assessed
#>     three candidates: stability not assessed
#>  the criterion does not fall along this line
#> 
#> 5 tables above, truncated to fit. get_results(x, what = ) returns any
#> of them whole, and get_results(x, what = "all") returns every one.
plot(models)
```
