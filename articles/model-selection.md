# Model selection: the evidence for the number of classes

In finite mixture modelling, the number of latent classes is not known a
priori and is not a parameter of the model; it must be determined by
comparing models fitted with different numbers of classes. No single
statistic identifies the correct number: information criteria,
likelihood-based tests and indices of classification quality each have
known limitations, and they frequently disagree (Nylund, Asparouhov and
Muthén, 2007). Class enumeration is therefore an inferential judgement
that integrates several sources of evidence and is reported
transparently. This guide describes each source of evidence that latents
computes for a set of candidate models: its definition, its
interpretation, and the conditions under which it is unreliable. It then
outlines how these sources are combined into a defensible and
reproducible decision. A complete applied example is given in [Class
enumeration in mixture
models](https://pak.dynasite.org/latents/articles/class-enumeration.md).

The evidence falls into five kinds, each answering a different question.

| Question | Evidence | Verbs |
|----|----|----|
| Which candidate balances fit and size best? | Information criteria and their minima | [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), `get_results(x, "best")`, `plot(x)` |
| Where do added classes stop paying off? | Changes between neighbouring candidates and the knee of each curve | `get_results(x, "deltas")`, `get_results(x, "knees")`, `plot(x, what = "deltas")`, `plot(x, what = "heatmap")` |
| Are the classes distinct? | Entropy, average posterior probability, odds of correct classification | `get_results(x, "classification")`, `plot(x, what = "diagnostics")` |
| Can the estimate be trusted? | Convergence, boundary estimates, replication of the best likelihood, the smallest class | `plot(x, what = "trust")`, `plot(x, what = "smallest")` |
| Is one more class supported? | The bootstrap likelihood ratio test of neighbouring candidates | [`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md) |

`plot(x, what = "selection")` draws the plots of the second, third and
fourth rows together, and
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
places chosen candidates side by side.

``` r

library(latents)
```

## The candidates

The examples use `course_engagement`: course enrolments nested in
students, with five standardized measures of activity per enrolment.
Every combination of one to six profiles and one to four group classes
is fitted with a covariance structure of free variances and zero
covariances (`"VVI"`). Three random starts per candidate keep the
example short; applied analyses use ten or more, because replication of
the best likelihood across starts is itself part of the evidence.

``` r

activity <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
models <- enumerate_classes(course_engagement, activity, "student",
                            n_profiles = 1:6, n_group_classes = 1:4,
                            model = "VVI", n_starts = 3, seed = 1)
models
#> Class enumeration: 24 candidate models
#>          candidate n_profiles n_group_classes model log_likelihood n_parameters   aic bic_groups
#>  profiles1_groups1          1               1   VVI         -10007           10 20035      20062
#>  profiles2_groups1          2               1   VVI          -8615           21 17272      17328
#>  profiles3_groups1          3               1   VVI          -8537           32 17137      17222
#>  profiles4_groups1          4               1   VVI          -8501           43 17088      17203
#>  profiles5_groups1          5               1   VVI          -8484           54 17077      17220
#>  profiles6_groups1          6               1   VVI          -8481           65 17093      17266
#>  profiles1_groups2          1               2   VVI             NA           NA    NA         NA
#>  profiles2_groups2          2               2   VVI          -8440           23 16926      16987
#>  profiles3_groups2          3               2   VVI          -8365           35 16799      16892
#>  profiles4_groups2          4               2   VVI          -8326           47 16746      16872
#>  profiles5_groups2          5               2   VVI          -8309           59 16735      16892
#>  profiles6_groups2          6               2   VVI          -8297           71 16737      16926
#>  profiles1_groups3          1               3   VVI             NA           NA    NA         NA
#>  profiles2_groups3          2               3   VVI          -8421           25 16892      16959
#>  profiles3_groups3          3               3   VVI          -8348           38 16773      16874
#>  profiles4_groups3          4               3   VVI          -8307           51 16716      16852
#>  profiles5_groups3          5               3   VVI          -8291           64 16709      16880
#>  profiles6_groups3          6               3   VVI          -8279           77 16712      16917
#>  profiles1_groups4          1               4   VVI             NA           NA    NA         NA
#>  profiles2_groups4          2               4   VVI          -8415           27 16884      16956
#>  profiles3_groups4          3               4   VVI          -8343           41 16767      16876
#>  profiles4_groups4          4               4   VVI          -8301           55 16711      16858
#>  profiles5_groups4          5               4   VVI          -8283           69 16704      16887
#>  profiles6_groups4          6               4   VVI          -8268           83 16703      16924
#>  bic_individual icl_individual profile_entropy group_entropy avepp_min avepp_max avepp_mean
#>           20088          20088              NA            NA        NA        NA         NA
#>           17383          17537           0.922            NA     0.974     0.982      0.978
#>           17306          18017           0.772            NA     0.832     0.964      0.886
#>           17314          18359           0.735            NA     0.803     0.849      0.835
#>           17361          18401           0.773            NA     0.804     0.962      0.861
#>           17435          19124           0.669            NA     0.627     0.812      0.729
#>              NA             NA              NA            NA        NA        NA         NA
#>           17047          17168           0.938         0.940     0.982     0.985      0.983
#>           16983          17640           0.790         0.943     0.849     0.979      0.893
#>           16994          18021           0.739         0.937     0.791     0.858      0.835
#>           17045          18076           0.775         0.965     0.803     0.860      0.839
#>           17110          18476           0.732         0.946     0.699     0.845      0.780
#>              NA             NA              NA            NA        NA        NA         NA
#>           17024          17142           0.940         0.816     0.982     0.986      0.984
#>           16973          17622           0.792         0.797     0.845     0.978      0.892
#>           16985          18014           0.739         0.823     0.801     0.854      0.836
#>           17046          18353           0.714         0.854     0.628     0.843      0.779
#>           17117          18780           0.674         0.853     0.655     0.795      0.743
#>              NA             NA              NA            NA        NA        NA         NA
#>           17026          17140           0.942         0.798     0.984     0.986      0.985
#>           16983          17630           0.793         0.772     0.843     0.980      0.892
#>           17000          18027           0.740         0.803     0.801     0.853      0.837
#>           17067          18207           0.751         0.830     0.786     0.853      0.820
#>           17139          18540           0.725         0.851     0.654     0.847      0.783
#>  avepp_below group_avepp_min group_avepp_max group_avepp_mean group_avepp_below converged boundary
#>           NA              NA              NA               NA                NA      TRUE    FALSE
#>            0              NA              NA               NA                NA      TRUE    FALSE
#>            0              NA              NA               NA                NA      TRUE    FALSE
#>            0              NA              NA               NA                NA      TRUE    FALSE
#>            0              NA              NA               NA                NA      TRUE    FALSE
#>            5              NA              NA               NA                NA      TRUE    FALSE
#>           NA              NA              NA               NA                NA     FALSE       NA
#>            0           0.969           0.992            0.980                 0      TRUE    FALSE
#>            0           0.977           0.992            0.985                 0      TRUE    FALSE
#>            1           0.975           0.990            0.983                 0      TRUE    FALSE
#>            0           0.990           0.996            0.993                 0      TRUE    FALSE
#>            4           0.972           0.987            0.979                 0      TRUE    FALSE
#>           NA              NA              NA               NA                NA     FALSE       NA
#>            0           0.900           0.962            0.923                 0      TRUE    FALSE
#>            0           0.864           0.977            0.913                 0      TRUE    FALSE
#>            0           0.874           0.981            0.923                 0      TRUE    FALSE
#>            1           0.917           0.987            0.941                 0      TRUE    FALSE
#>            6           0.915           0.984            0.941                 0      TRUE    FALSE
#>           NA              NA              NA               NA                NA     FALSE       NA
#>            0           0.851           0.912            0.877                 0      TRUE    FALSE
#>            0           0.841           0.902            0.864                 0      TRUE    FALSE
#>            0           0.867           0.903            0.885                 0      TRUE    FALSE
#>            1           0.882           0.926            0.904                 0      TRUE    FALSE
#>            3           0.898           0.940            0.916                 0      TRUE    FALSE
#>  n_best_replicated
#>                  3
#>                  3
#>                  3
#>                  3
#>                  1
#>                  1
#>                 NA
#>                  3
#>                  3
#>                  3
#>                  2
#>                  1
#>                 NA
#>                  3
#>                  3
#>                  3
#>                  2
#>                  2
#>                 NA
#>                  3
#>                  3
#>                  2
#>                  1
#>                  1
#> 3 candidate(s) did not converge; see as.data.frame(x).
#> No candidate is selected automatically. Compare one convention consistently.
```

The grid holds 24 candidates, of which 21 were fitted. A candidate with
one profile and several group classes is not identified, because group
classes that share a single profile cannot differ, so it is listed and
left unfitted. Every candidate is named by its numbers of classes,
`profiles4_groups3` for four profiles and three group classes, and the
verbs below accept that name.

## Before fitting: starts, scaling and missing values

Three choices are made before any candidate is fitted. Each is passed
through the `...` of
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md),
[`enumerate_lpa()`](https://pak.dynasite.org/latents/reference/enumerate_lpa.md)
and
[`enumerate_lca()`](https://pak.dynasite.org/latents/reference/enumerate_lca.md),
so every candidate receives the same setting, and each is reported with
the choice of the number of classes.

### Random starts

`n_starts` sets the number of random starts per candidate. More starts
make it less likely that a candidate is reported at a local maximum of
its likelihood, and a candidate stuck at a local maximum usually
separates its classes worse than it would at the global one: refitting
with more starts often raises its entropy and AvePP. The whole grid is
refitted with the same, larger number of starts, not only the candidates
of interest; otherwise the refitted candidates gain likelihood and
separation that the others were not given the chance to reach.

``` r

models <- enumerate_classes(course_engagement, activity, "student",
                            n_profiles = 1:6, n_group_classes = 1:4,
                            model = "VVI", n_starts = 20, seed = 1)
```

### Centring and scaling

`centering` subtracts a mean and `scaling` divides by a standard
deviation, each at one of two levels: `"grand"` uses one value per
indicator over all observations, and `"person"` one per group. Together
they give the classical z-scores:
`centering = "grand", scaling = "grand"` standardizes each indicator
over the sample, and `centering = "person", scaling = "person"`
standardizes it within each group, so the profiles describe how an
observation departs from its own group’s usual level and spread.

For the covariance structures that give every indicator its own variance
(`EEI`, `VVI`, `EEE`, `VVV`, `EVI`, `VEI`), grand scaling changes only
the units of the estimates: memberships, entropy and AvePP are
unchanged. The spherical `EII` and `VII` force one variance on every
indicator, so an unscaled indicator with a wide range would dominate the
classes; with them, the indicators are scaled first.

``` r

srl_indicators <- c("cognitive_strategies", "intrinsic_value",
                    "self_efficacy", "self_regulation", "test_anxiety")
spherical <- enumerate_lpa(srl, srl_indicators, n_profiles = 2:4,
                           model = "EII", centering = "grand",
                           scaling = "grand", n_starts = 5, seed = 1)
get_results(spherical, "best")
#>    criterion  convention candidate n_profiles n_group_classes model value
#> 1        aic        <NA> profiles4          4               1   EII  3372
#> 2        kic        <NA> profiles4          4               1   EII  3399
#> 3        bic      groups profiles4          4               1   EII  3461
#> 4        bic individuals profiles4          4               1   EII  3461
#> 5      sabic      groups profiles4          4               1   EII  3385
#> 6      sabic individuals profiles4          4               1   EII  3385
#> 7       caic      groups profiles4          4               1   EII  3485
#> 8       caic individuals profiles4          4               1   EII  3485
#> 9        awe      groups profiles4          4               1   EII  3670
#> 10       awe individuals profiles4          4               1   EII  3779
#> 11       icl      groups profiles4          4               1   EII  3461
#> 12       icl individuals profiles4          4               1   EII  3570
#> 13       clc      groups profiles4          4               1   EII  3324
#> 14       clc individuals profiles4          4               1   EII  3433
```

Scaling moves the log likelihood, so information criteria compare only
candidates fitted with the same centring and scaling, as the candidates
of one grid are.

### Missing values

With `missing = "fiml"`, a row with some indicators missing is fitted
from the indicators it has, and nothing is imputed. A row with every
indicator missing adds nothing to the likelihood, and its posterior is
only the prior its group implies: it was never classified. Such rows are
left out of the sample size of the individual-level criteria and out of
entropy, AvePP and the classification tables, as Mplus leaves them out
of its analysis; the fit’s `n_informative` counts the rows that remain,
and its print reports how many were left out. A group whose every row
has nothing observed is left out of the group sample size and the
group-level statistics in the same way (`n_informative_groups`); a group
with at least one observed row stays, and the effective class counts are
of the classified rows and groups. Rows with some values missing stay in
every statistic: they are classified, less certainly, from fewer
indicators. `drop_empty = FALSE` restores the counting of latents 0.9.37
and earlier, so that results reported with those versions can be
reproduced.

``` r

set.seed(3)
incomplete <- course_engagement
empty <- sample(nrow(incomplete), 150)
incomplete[empty, activity] <- NA
partial <- setdiff(sample(nrow(incomplete), 150), empty)
incomplete[partial, "attendance"] <- NA
incomplete_fit <- multilpa(incomplete, activity, "student", n_profiles = 3,
                           n_group_classes = 2, model = "VVI",
                           missing = "fiml", n_starts = 3, seed = 1)
get_results(incomplete_fit, "entropy")
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         3    1272      321.52            0.770
#> 2      groups         2     106        5.42            0.926
```

The entropy of the individuals is taken over 1272 rows: the 1,422
enrolments less the 150 with nothing observed. Dropping those rows from
the data by hand gives the same estimates, likelihood and classification
statistics; dropping rows with only some values missing would change the
sample and, if their missingness depends on what was observed, bias the
profiles.

## Information criteria

### Definitions

Each criterion starts from the deviance, `-2L`, where `L` is the
maximized log likelihood, and adds a charge. With `q` free parameters,
`n` the sample size and `EN` the classification entropy,
`EN = -sum(p log p)` summed over every unit and class of the posterior
probabilities `p`:

| Criterion | Formula | Charge depends on | Reference |
|----|----|----|----|
| `aic` | `-2L + 2q` | parameters only | Akaike (1974) |
| `kic` | `-2L + 3(q + 1)` | parameters only | Cavanaugh (1999) |
| `bic` | `-2L + q log(n)` | sample size | Schwarz (1978) |
| `sabic` | `-2L + q log((n + 2) / 24)` | sample size | Sclove (1987) |
| `caic` | `-2L + q (log(n) + 1)` | sample size | Bozdogan (1987) |
| `icl` | `-2L + q log(n) + 2 EN` | sample size and classification | Biernacki, Celeux and Govaert (2000) |
| `awe` | `-2(L - EN) + 2q (1.5 + log(n))` | sample size and classification | Banfield and Raftery (1993) |
| `clc` | `-2L + 2 EN` | classification only | Biernacki and Govaert (1997) |

Lower is better for every criterion. The charge per parameter is 2 for
AIC, 3 for KIC, `log((n + 2) / 24)` for SABIC, `log(n)` for BIC and
`log(n) + 1` for CAIC, so above a few dozen units AIC is the most
lenient and CAIC the strictest of these, and a larger charge favours
fewer classes. ICL and AWE add a charge for classification uncertainty,
because a class that overlaps its neighbours raises `EN`, and so favour
fewer and better separated classes than BIC. CLC has no charge per
parameter, only the entropy charge, so it can favour many classes when
they are well separated and few when they are not. CLC uses the entropy
sum as ICL and AWE do. Some software reports a CLC built from relative
entropy, which is bounded by one and charges almost nothing, and the two
numbers are not comparable.

### Two conventions for the sample size

In a two-level model the sample size is either the number of groups,
here the 106 students, or the number of observations, the 1,422
enrolments. latents reports every criterion that depends on a sample
size under both conventions: the columns ending `_groups` count groups
and use the posterior probabilities of the group classes for `EN`, and
the columns ending `_individual` count observations and use the
posterior probabilities of the profiles. Neither is the definition of a
multilevel BIC; they are two stated choices, and an analysis picks one
and keeps it for every candidate. AIC and KIC have no sample size and
carry no convention.

In a single-level grid there are no groups, the two sample sizes
coincide, and BIC, SABIC and CAIC are the same under both conventions.
The entropy criteria are not: their `_groups` columns use the
classification of the group classes, of which a single-level model has
one, so they carry no classification charge and ICL there equals BIC. In
a single-level grid, read the `_individual` columns of ICL, AWE and CLC.

### The table of every candidate

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per candidate with every criterion under both conventions,
beside its log likelihood, number of parameters, classification
summaries and estimation checks.

``` r

as.data.frame(models)
#>            candidate n_profiles n_group_classes model log_likelihood n_parameters   aic   kic
#> 1  profiles1_groups1          1               1   VVI         -10007           10 20035 20048
#> 2  profiles2_groups1          2               1   VVI          -8615           21 17272 17296
#> 3  profiles3_groups1          3               1   VVI          -8537           32 17137 17172
#> 4  profiles4_groups1          4               1   VVI          -8501           43 17088 17134
#> 5  profiles5_groups1          5               1   VVI          -8484           54 17077 17134
#> 6  profiles6_groups1          6               1   VVI          -8481           65 17093 17161
#> 7  profiles1_groups2          1               2   VVI             NA           NA    NA    NA
#> 8  profiles2_groups2          2               2   VVI          -8440           23 16926 16952
#> 9  profiles3_groups2          3               2   VVI          -8365           35 16799 16837
#> 10 profiles4_groups2          4               2   VVI          -8326           47 16746 16796
#> 11 profiles5_groups2          5               2   VVI          -8309           59 16735 16797
#> 12 profiles6_groups2          6               2   VVI          -8297           71 16737 16811
#> 13 profiles1_groups3          1               3   VVI             NA           NA    NA    NA
#> 14 profiles2_groups3          2               3   VVI          -8421           25 16892 16920
#> 15 profiles3_groups3          3               3   VVI          -8348           38 16773 16814
#> 16 profiles4_groups3          4               3   VVI          -8307           51 16716 16770
#> 17 profiles5_groups3          5               3   VVI          -8291           64 16709 16776
#> 18 profiles6_groups3          6               3   VVI          -8279           77 16712 16792
#> 19 profiles1_groups4          1               4   VVI             NA           NA    NA    NA
#> 20 profiles2_groups4          2               4   VVI          -8415           27 16884 16914
#> 21 profiles3_groups4          3               4   VVI          -8343           41 16767 16811
#> 22 profiles4_groups4          4               4   VVI          -8301           55 16711 16769
#> 23 profiles5_groups4          5               4   VVI          -8283           69 16704 16776
#> 24 profiles6_groups4          6               4   VVI          -8268           83 16703 16789
#>    bic_groups bic_individual sabic_groups sabic_individual caic_groups caic_individual awe_groups
#> 1       20062          20088        20030            20056       20072           20098      20138
#> 2       17328          17383        17262            17316       17349           17404      17489
#> 3       17222          17306        17121            17204       17254           17338      17468
#> 4       17203          17314        17067            17178       17246           17357      17532
#> 5       17220          17361        17050            17189       17274           17415      17634
#> 6       17266          17435        17061            17228       17331           17500      17764
#> 7          NA             NA           NA               NA          NA              NA         NA
#> 8       16987          17047        16914            16974       17010           17070      17172
#> 9       16892          16983        16782            16872       16927           17018      17169
#> 10      16872          16994        16723            16844       16919           17041      17241
#> 11      16892          17045        16706            16858       16951           17104      17350
#> 12      16926          17110        16701            16884       16997           17181      17478
#> 13         NA             NA           NA               NA          NA              NA         NA
#> 14      16959          17024        16880            16944       16984           17049      17193
#> 15      16874          16973        16754            16852       16912           17011      17212
#> 16      16852          16985        16691            16823       16903           17036      17284
#> 17      16880          17046        16677            16842       16944           17110      17404
#> 18      16917          17117        16674            16872       16994           17194      17541
#> 19         NA             NA           NA               NA          NA              NA         NA
#> 20      16956          17026        16871            16941       16983           17053      17223
#> 21      16876          16983        16747            16853       16917           17024      17258
#> 22      16858          17000        16684            16826       16913           17055      17337
#> 23      16887          17067        16669            16847       16956           17136      17466
#> 24      16924          17139        16661            16876       17007           17222      17604
#>    awe_individual icl_groups icl_individual clc_groups clc_individual profile_entropy group_entropy
#> 1           20190      20062          20088      20015          20015              NA            NA
#> 2           17753      17328          17537      17230          17385           0.922            NA
#> 3           18345      17222          18017      17073          17785           0.772            NA
#> 4           18800      17203          18359      17002          18047           0.735            NA
#> 5           18955      17220          18401      16969          18009           0.773            NA
#> 6           19791      17266          19124      16963          18652           0.669            NA
#> 7              NA         NA             NA         NA             NA              NA            NA
#> 8           17404      16996          17168      16889          17002           0.938         0.940
#> 9           18000      16901          17640      16737          17386           0.790         0.943
#> 10          18503      16881          18021      16662          17680           0.739         0.937
#> 11          18682      16897          18076      16622          17648           0.775         0.965
#> 12          19205      16934          18476      16602          17961           0.732         0.946
#> 13             NA         NA             NA         NA             NA              NA            NA
#> 14          17398      17001          17142      16885          16960           0.940         0.816
#> 15          18012      16921          17622      16744          17346           0.792         0.797
#> 16          18537      16893          18014      16656          17644           0.739         0.823
#> 17          19010      16914          18353      16615          17888           0.714         0.854
#> 18          19570      16951          18780      16592          18221           0.674         0.853
#> 19             NA         NA             NA         NA             NA              NA            NA
#> 20          17417      17016          17140      16890          16944           0.942         0.798
#> 21          18051      16943          17630      16752          17332           0.793         0.772
#> 22          18592      16916          18027      16659          17628           0.740         0.803
#> 23          18915      16937          18207      16616          17706           0.751         0.830
#> 24          19392      16967          18540      16580          17938           0.725         0.851
#>    avepp_min avepp_max avepp_mean avepp_below group_avepp_min group_avepp_max group_avepp_mean
#> 1         NA        NA         NA          NA              NA              NA               NA
#> 2      0.974     0.982      0.978           0              NA              NA               NA
#> 3      0.832     0.964      0.886           0              NA              NA               NA
#> 4      0.803     0.849      0.835           0              NA              NA               NA
#> 5      0.804     0.962      0.861           0              NA              NA               NA
#> 6      0.627     0.812      0.729           5              NA              NA               NA
#> 7         NA        NA         NA          NA              NA              NA               NA
#> 8      0.982     0.985      0.983           0           0.969           0.992            0.980
#> 9      0.849     0.979      0.893           0           0.977           0.992            0.985
#> 10     0.791     0.858      0.835           1           0.975           0.990            0.983
#> 11     0.803     0.860      0.839           0           0.990           0.996            0.993
#> 12     0.699     0.845      0.780           4           0.972           0.987            0.979
#> 13        NA        NA         NA          NA              NA              NA               NA
#> 14     0.982     0.986      0.984           0           0.900           0.962            0.923
#> 15     0.845     0.978      0.892           0           0.864           0.977            0.913
#> 16     0.801     0.854      0.836           0           0.874           0.981            0.923
#> 17     0.628     0.843      0.779           1           0.917           0.987            0.941
#> 18     0.655     0.795      0.743           6           0.915           0.984            0.941
#> 19        NA        NA         NA          NA              NA              NA               NA
#> 20     0.984     0.986      0.985           0           0.851           0.912            0.877
#> 21     0.843     0.980      0.892           0           0.841           0.902            0.864
#> 22     0.801     0.853      0.837           0           0.867           0.903            0.885
#> 23     0.786     0.853      0.820           1           0.882           0.926            0.904
#> 24     0.654     0.847      0.783           3           0.898           0.940            0.916
#>    group_avepp_below converged boundary n_best_replicated warnings
#> 1                 NA      TRUE    FALSE                 3         
#> 2                 NA      TRUE    FALSE                 3         
#> 3                 NA      TRUE    FALSE                 3         
#> 4                 NA      TRUE    FALSE                 3         
#> 5                 NA      TRUE    FALSE                 1         
#> 6                 NA      TRUE    FALSE                 1         
#> 7                 NA     FALSE       NA                NA         
#> 8                  0      TRUE    FALSE                 3         
#> 9                  0      TRUE    FALSE                 3         
#> 10                 0      TRUE    FALSE                 3         
#> 11                 0      TRUE    FALSE                 2         
#> 12                 0      TRUE    FALSE                 1         
#> 13                NA     FALSE       NA                NA         
#> 14                 0      TRUE    FALSE                 3         
#> 15                 0      TRUE    FALSE                 3         
#> 16                 0      TRUE    FALSE                 3         
#> 17                 0      TRUE    FALSE                 2         
#> 18                 0      TRUE    FALSE                 2         
#> 19                NA     FALSE       NA                NA         
#> 20                 0      TRUE    FALSE                 3         
#> 21                 0      TRUE    FALSE                 3         
#> 22                 0      TRUE    FALSE                 2         
#> 23                 0      TRUE    FALSE                 1         
#> 24                 0      TRUE    FALSE                 1         
#>                                                                                     error
#> 1                                                                                    <NA>
#> 2                                                                                    <NA>
#> 3                                                                                    <NA>
#> 4                                                                                    <NA>
#> 5                                                                                    <NA>
#> 6                                                                                    <NA>
#> 7  Multiple group classes are not identifiable with one profile or only singleton groups.
#> 8                                                                                    <NA>
#> 9                                                                                    <NA>
#> 10                                                                                   <NA>
#> 11                                                                                   <NA>
#> 12                                                                                   <NA>
#> 13 Multiple group classes are not identifiable with one profile or only singleton groups.
#> 14                                                                                   <NA>
#> 15                                                                                   <NA>
#> 16                                                                                   <NA>
#> 17                                                                                   <NA>
#> 18                                                                                   <NA>
#> 19 Multiple group classes are not identifiable with one profile or only singleton groups.
#> 20                                                                                   <NA>
#> 21                                                                                   <NA>
#> 22                                                                                   <NA>
#> 23                                                                                   <NA>
#> 24                                                                                   <NA>
```

### The minimum of each criterion

`get_results(x, "best")` reduces that table to one row per criterion and
convention: the candidate at the criterion’s minimum, its numbers of
classes and the minimum value. It answers “which candidate does each
criterion favour?” and nothing more.

``` r

get_results(models, "best")
#>    criterion  convention         candidate n_profiles n_group_classes model value
#> 1        aic        <NA> profiles6_groups4          6               4   VVI 16703
#> 2        kic        <NA> profiles4_groups4          4               4   VVI 16769
#> 3        bic      groups profiles4_groups3          4               3   VVI 16852
#> 4        bic individuals profiles3_groups3          3               3   VVI 16973
#> 5      sabic      groups profiles6_groups4          6               4   VVI 16661
#> 6      sabic individuals profiles4_groups3          4               3   VVI 16823
#> 7       caic      groups profiles4_groups3          4               3   VVI 16903
#> 8       caic individuals profiles3_groups3          3               3   VVI 17011
#> 9        awe      groups profiles3_groups2          3               2   VVI 17169
#> 10       awe individuals profiles2_groups3          2               3   VVI 17398
#> 11       icl      groups profiles4_groups2          4               2   VVI 16881
#> 12       icl individuals profiles2_groups4          2               4   VVI 17140
#> 13       clc      groups profiles6_groups4          6               4   VVI 16580
#> 14       clc individuals profiles2_groups4          2               4   VVI 16944
```

Here the criteria name 8 different candidates. BIC with students as the
units favours profiles4_groups3, BIC with enrolments profiles3_groups3,
ICL with students profiles4_groups2, and AIC profiles6_groups4.
Disagreement is the usual outcome, not a failure: the criteria were
built for different purposes, and the size of their charges explains
much of the disagreement: AIC, with the lightest charge, names the
largest candidate, and the entropy criteria name smaller ones.

A minimum on the edge of the grid, at the largest number of profiles or
group classes fitted, means the criterion was still falling there. Such
a minimum says only that the range was too narrow. The range is widened
until the minimum lies inside it, or, when the criterion keeps falling
at every size the data can support, the size of each decrease becomes
the evidence, as in the next section.

### Plotting the criteria

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
criteria against the number of profiles, one line per number of group
classes, and rings the lowest value of each. `what` names criteria by
their columns, and `what = "all"` draws every one.

``` r

plot(models)
```

![](model-selection_files/figure-html/plot-criteria-1.png)

``` r

plot(models, what = "bic_groups")
```

![](model-selection_files/figure-html/plot-bic-1.png)

## Where added classes stop paying off

### Deltas

When a criterion falls across the whole range, its minimum carries
little information, but the size of each fall does.
`get_results(x, "deltas")` reports, for every candidate, the decrease in
a criterion from its neighbour with one profile fewer (`delta_profile`,
from `from_profile`) and from its neighbour with one group class fewer
(`delta_group_class`, from `from_group_class`). A positive delta means
the added class lowered the criterion; a negative delta means it raised
it. A candidate with no such neighbour in the grid has `NA`. The
criterion is BIC with groups as the units unless `using` names another
column.

``` r

get_results(models, "deltas")
#>            candidate n_profiles n_group_classes model  criterion value delta_profile
#> 1  profiles1_groups1          1               1   VVI bic_groups 20062            NA
#> 2  profiles2_groups1          2               1   VVI bic_groups 17328        2733.2
#> 3  profiles3_groups1          3               1   VVI bic_groups 17222         105.9
#> 4  profiles4_groups1          4               1   VVI bic_groups 17203          19.9
#> 5  profiles5_groups1          5               1   VVI bic_groups 17220         -17.9
#> 6  profiles6_groups1          6               1   VVI bic_groups 17266         -45.6
#> 7  profiles1_groups2          1               2   VVI bic_groups    NA            NA
#> 8  profiles2_groups2          2               2   VVI bic_groups 16987            NA
#> 9  profiles3_groups2          3               2   VVI bic_groups 16892          94.6
#> 10 profiles4_groups2          4               2   VVI bic_groups 16872          20.7
#> 11 profiles5_groups2          5               2   VVI bic_groups 16892         -20.6
#> 12 profiles6_groups2          6               2   VVI bic_groups 16926         -33.5
#> 13 profiles1_groups3          1               3   VVI bic_groups    NA            NA
#> 14 profiles2_groups3          2               3   VVI bic_groups 16959            NA
#> 15 profiles3_groups3          3               3   VVI bic_groups 16874          84.6
#> 16 profiles4_groups3          4               3   VVI bic_groups 16852          21.9
#> 17 profiles5_groups3          5               3   VVI bic_groups 16880         -27.5
#> 18 profiles6_groups3          6               3   VVI bic_groups 16917         -37.5
#> 19 profiles1_groups4          1               4   VVI bic_groups    NA            NA
#> 20 profiles2_groups4          2               4   VVI bic_groups 16956            NA
#> 21 profiles3_groups4          3               4   VVI bic_groups 16876          79.9
#> 22 profiles4_groups4          4               4   VVI bic_groups 16858          18.8
#> 23 profiles5_groups4          5               4   VVI bic_groups 16887         -29.7
#> 24 profiles6_groups4          6               4   VVI bic_groups 16924         -36.3
#>         from_profile delta_group_class  from_group_class
#> 1               <NA>                NA              <NA>
#> 2  profiles1_groups1                NA              <NA>
#> 3  profiles2_groups1                NA              <NA>
#> 4  profiles3_groups1                NA              <NA>
#> 5  profiles4_groups1                NA              <NA>
#> 6  profiles5_groups1                NA              <NA>
#> 7               <NA>                NA              <NA>
#> 8               <NA>            341.51 profiles2_groups1
#> 9  profiles2_groups2            330.20 profiles3_groups1
#> 10 profiles3_groups2            330.96 profiles4_groups1
#> 11 profiles4_groups2            328.23 profiles5_groups1
#> 12 profiles5_groups2            340.42 profiles6_groups1
#> 13              <NA>                NA              <NA>
#> 14              <NA>             28.18 profiles2_groups2
#> 15 profiles2_groups3             18.19 profiles3_groups2
#> 16 profiles3_groups3             19.43 profiles4_groups2
#> 17 profiles4_groups3             12.55 profiles5_groups2
#> 18 profiles5_groups3              8.51 profiles6_groups2
#> 19              <NA>                NA              <NA>
#> 20              <NA>              2.37 profiles2_groups3
#> 21 profiles2_groups4             -2.43 profiles3_groups3
#> 22 profiles3_groups4             -5.54 profiles4_groups3
#> 23 profiles4_groups4             -7.79 profiles5_groups3
#> 24 profiles5_groups4             -6.62 profiles6_groups3
```

With three group classes, the third profile lowers BIC by 85 and the
fourth by 22, and the fifth raises it by 27, which the table shows as a
negative decrease. The candidate with two profiles has no delta along
the profiles, because its neighbour with one profile and three group
classes is not identified. Shrinking falls mark the point after which
further profiles add detail rather than structure, and the first
negative delta marks where they stop improving the criterion at all. The
deltas are plotted for both directions, along the number of profiles on
the left and along the number of group classes on the right.

``` r

plot(models, what = "deltas")
```

![](model-selection_files/figure-html/plot-deltas-1.png)

Any criterion can be read this way:

``` r

plot(models, what = "deltas", using = "icl_individual")
```

![](model-selection_files/figure-html/plot-deltas-icl-1.png)

### Knees

The knee, or elbow, is where a falling curve bends from steep to flat.
`get_results(x, "knees")` locates it by a rule inspired by Kneedle
(Satopää, Albrecht, Irwin and Raghavan, 2011), simplified to its
geometric core:

1.  both axes of the curve, the number of classes and the criterion, are
    rescaled to the unit interval;
2.  a straight line is drawn from the first candidate to the last;
3.  the knee is the candidate furthest below that line, measured
    perpendicular to it, and `distance` is that perpendicular distance
    in the rescaled units.

Kneedle’s smoothing and sensitivity thresholds are left out, so a nearly
straight curve can still name a knee with a very small `distance`. The
distance therefore reports the strength of the bend as well as its
place. A curve whose last value is not below its first has no knee,
because a bend in a rising criterion marks no diminishing return.

Knees are found along three kinds of line: along the number of profiles
at each number of group classes (`direction = "profiles"`), along the
number of group classes at each number of profiles (`"group_classes"`),
and across all candidates ordered by their number of parameters
(`"parameters"`). The columns `from` and `to` give the ends of the line
and `n_points` its length.

``` r

get_results(models, "knees")
#>        direction            line  criterion              from                to n_points
#> 1       profiles   1 group class bic_groups profiles1_groups1 profiles6_groups1        6
#> 2       profiles 2 group classes bic_groups profiles2_groups2 profiles6_groups2        5
#> 3       profiles 3 group classes bic_groups profiles2_groups3 profiles6_groups3        5
#> 4       profiles 4 group classes bic_groups profiles2_groups4 profiles6_groups4        5
#> 5  group_classes       1 profile bic_groups profiles1_groups1 profiles1_groups1        1
#> 6  group_classes      2 profiles bic_groups profiles2_groups1 profiles2_groups4        4
#> 7  group_classes      3 profiles bic_groups profiles3_groups1 profiles3_groups4        4
#> 8  group_classes      4 profiles bic_groups profiles4_groups1 profiles4_groups4        4
#> 9  group_classes      5 profiles bic_groups profiles5_groups1 profiles5_groups4        4
#> 10 group_classes      6 profiles bic_groups profiles6_groups1 profiles6_groups4        4
#> 11    parameters  all candidates bic_groups profiles1_groups1 profiles6_groups4       21
#>       knee_candidate distance knee_without_largest stable                               note
#> 1  profiles2_groups1    0.544    profiles2_groups1   TRUE                                   
#> 2  profiles4_groups2    0.648    profiles3_groups2  FALSE moves when the largest is left out
#> 3  profiles4_groups3    0.750    profiles3_groups3  FALSE moves when the largest is left out
#> 4  profiles4_groups4    0.792    profiles3_groups4  FALSE moves when the largest is left out
#> 5               <NA>       NA                 <NA>     NA        fewer than three candidates
#> 6  profiles2_groups2    0.413    profiles2_groups2   TRUE                                   
#> 7  profiles3_groups2    0.438    profiles3_groups2   TRUE                                   
#> 8  profiles4_groups2    0.439    profiles4_groups2   TRUE                                   
#> 9  profiles5_groups2    0.456    profiles5_groups2   TRUE                                   
#> 10 profiles6_groups2    0.463    profiles6_groups2   TRUE                                   
#> 11 profiles2_groups2    0.561    profiles2_groups2   TRUE
```

The straight line ends at the largest candidate fitted, so the knee
depends on where the grid stops: fitting one more candidate moves the
line and can move the knee. The table tests that directly.
`knee_without_largest` is the knee of the same line with its last
candidate removed, and `stable` is `TRUE` when the two agree. A knee
that moves when the range shrinks by one candidate is weak evidence.
`note` explains every row without a full answer:

| `note` | Meaning |
|----|----|
| empty | a knee was found and is stable |
| `fewer than three candidates` | a line needs three points for a knee |
| `three candidates: stability not assessed` | `stable` needs four points |
| `the criterion does not fall along this line` | no diminishing return to find |
| `the curve does not bend below the line` | the curve is straight or bends upward |
| `no knee when the largest is left out` | the shortened line has no knee |
| `moves when the largest is left out` | the knee depends on the range fitted |

With three group classes, the knee along the number of profiles is
profiles4_groups3 at a distance of 0.75. Without the largest candidate
it is profiles3_groups3, so `stable` is FALSE: the bend lies between
three and four profiles, and which of the two the rule names depends on
where the grid stops. The deltas locate the same region with less
ambiguity, which is why the two are read together.

### The criterion over the whole grid

The heatmap shows every candidate at once, coloured by how far its
criterion lies above the lowest value in the grid. It reads the same
criterion as the deltas and takes the same `using`.

``` r

plot(models, what = "heatmap")
```

![](model-selection_files/figure-html/heatmap-1.png)

## Are the classes distinct?

A candidate can fit well and still separate its classes poorly, leaving
many units between classes and their assignments uncertain. Three
statistics describe the separation, each from the posterior
probabilities of membership.

- **Relative entropy**, `1 - EN / (n log K)` for `K` classes, runs from
  0, every unit equally likely to belong to every class, to 1, every
  unit certain. It summarizes a whole level in one number and is `NA`
  for a level with one class, where it is undefined rather than perfect.
- **Average posterior probability (AvePP)** of a class is the mean
  posterior probability of that class among the units assigned to it.
  Values of 0.7 to 0.8 and above are commonly read as adequate
  separation (Nagin, 2005). It is reported per class, so a single poorly
  separated class is visible where entropy would average it away.
- **Odds of correct classification (OCC)** of a class is
  `(p / (1 - p)) / (r / (1 - r))`, with `p` its AvePP and `r` its
  estimated share. It compares classification against a guess from the
  class size alone: near 1 is no better than the guess, and 5 or more is
  the usual floor for adequate separation (Nagin, 2005).

The table of candidates summarizes each level by its relative entropy
and the smallest, largest and mean AvePP, and counts the classes below a
cutoff of 0.8, set with `cutoff`. `get_results(x, "classification")`
gives every class of every candidate: its size by modal assignment and
by estimate, its AvePP, its OCC and the relative entropy of its level.

``` r

get_results(models, "classification")
#>             candidate n_profiles n_group_classes model       level class n_modal proportion_modal
#> 1   profiles2_groups1          2               1   VVI individuals     1     833          0.58579
#> 2   profiles2_groups1          2               1   VVI individuals     2     589          0.41421
#> 3   profiles3_groups1          3               1   VVI individuals     1     452          0.31786
#> 4   profiles3_groups1          3               1   VVI individuals     2     555          0.39030
#> 5   profiles3_groups1          3               1   VVI individuals     3     415          0.29184
#> 6   profiles4_groups1          4               1   VVI individuals     1     481          0.33826
#> 7   profiles4_groups1          4               1   VVI individuals     2     157          0.11041
#> 8   profiles4_groups1          4               1   VVI individuals     3     419          0.29466
#> 9   profiles4_groups1          4               1   VVI individuals     4     365          0.25668
#> 10  profiles5_groups1          5               1   VVI individuals     1     355          0.24965
#> 11  profiles5_groups1          5               1   VVI individuals     2     419          0.29466
#> 12  profiles5_groups1          5               1   VVI individuals     3     486          0.34177
#> 13  profiles5_groups1          5               1   VVI individuals     4       5          0.00352
#> 14  profiles5_groups1          5               1   VVI individuals     5     157          0.11041
#> 15  profiles6_groups1          6               1   VVI individuals     1     372          0.26160
#> 16  profiles6_groups1          6               1   VVI individuals     2     157          0.11041
#> 17  profiles6_groups1          6               1   VVI individuals     3     378          0.26582
#> 18  profiles6_groups1          6               1   VVI individuals     4      62          0.04360
#> 19  profiles6_groups1          6               1   VVI individuals     5     337          0.23699
#> 20  profiles6_groups1          6               1   VVI individuals     6     116          0.08158
#> 21  profiles2_groups2          2               2   VVI individuals     1     836          0.58790
#> 22  profiles2_groups2          2               2   VVI individuals     2     586          0.41210
#> 23  profiles2_groups2          2               2   VVI      groups     1      71          0.66981
#> 24  profiles2_groups2          2               2   VVI      groups     2      35          0.33019
#> 25  profiles3_groups2          3               2   VVI individuals     1     472          0.33193
#> 26  profiles3_groups2          3               2   VVI individuals     2     570          0.40084
#> 27  profiles3_groups2          3               2   VVI individuals     3     380          0.26723
#> 28  profiles3_groups2          3               2   VVI      groups     1      71          0.66981
#> 29  profiles3_groups2          3               2   VVI      groups     2      35          0.33019
#> 30  profiles4_groups2          4               2   VVI individuals     1     490          0.34459
#> 31  profiles4_groups2          4               2   VVI individuals     2     166          0.11674
#> 32  profiles4_groups2          4               2   VVI individuals     3     418          0.29395
#> 33  profiles4_groups2          4               2   VVI individuals     4     348          0.24473
#> 34  profiles4_groups2          4               2   VVI      groups     1      71          0.66981
#> 35  profiles4_groups2          4               2   VVI      groups     2      35          0.33019
#> 36  profiles5_groups2          5               2   VVI individuals     1     483          0.33966
#> 37  profiles5_groups2          5               2   VVI individuals     2     423          0.29747
#> 38  profiles5_groups2          5               2   VVI individuals     3      31          0.02180
#> 39  profiles5_groups2          5               2   VVI individuals     4     327          0.22996
#> 40  profiles5_groups2          5               2   VVI individuals     5     158          0.11111
#> 41  profiles5_groups2          5               2   VVI      groups     1      72          0.67925
#> 42  profiles5_groups2          5               2   VVI      groups     2      34          0.32075
#> 43  profiles6_groups2          6               2   VVI individuals     1     465          0.32700
#> 44  profiles6_groups2          6               2   VVI individuals     2     171          0.12025
#> 45  profiles6_groups2          6               2   VVI individuals     3     317          0.22293
#> 46  profiles6_groups2          6               2   VVI individuals     4      90          0.06329
#> 47  profiles6_groups2          6               2   VVI individuals     5     343          0.24121
#> 48  profiles6_groups2          6               2   VVI individuals     6      36          0.02532
#> 49  profiles6_groups2          6               2   VVI      groups     1      69          0.65094
#> 50  profiles6_groups2          6               2   VVI      groups     2      37          0.34906
#> 51  profiles2_groups3          2               3   VVI individuals     1     832          0.58509
#> 52  profiles2_groups3          2               3   VVI individuals     2     590          0.41491
#> 53  profiles2_groups3          2               3   VVI      groups     1      25          0.23585
#> 54  profiles2_groups3          2               3   VVI      groups     2      31          0.29245
#> 55  profiles2_groups3          2               3   VVI      groups     3      50          0.47170
#> 56  profiles3_groups3          3               3   VVI individuals     1     367          0.25809
#> 57  profiles3_groups3          3               3   VVI individuals     2     576          0.40506
#> 58  profiles3_groups3          3               3   VVI individuals     3     479          0.33685
#> 59  profiles3_groups3          3               3   VVI      groups     1      25          0.23585
#> 60  profiles3_groups3          3               3   VVI      groups     2      30          0.28302
#> 61  profiles3_groups3          3               3   VVI      groups     3      51          0.48113
#> 62  profiles4_groups3          4               3   VVI individuals     1     508          0.35724
#> 63  profiles4_groups3          4               3   VVI individuals     2     402          0.28270
#> 64  profiles4_groups3          4               3   VVI individuals     3     192          0.13502
#> 65  profiles4_groups3          4               3   VVI individuals     4     320          0.22504
#> 66  profiles4_groups3          4               3   VVI      groups     1      28          0.26415
#> 67  profiles4_groups3          4               3   VVI      groups     2      48          0.45283
#> 68  profiles4_groups3          4               3   VVI      groups     3      30          0.28302
#> 69  profiles5_groups3          5               3   VVI individuals     1     123          0.08650
#> 70  profiles5_groups3          5               3   VVI individuals     2     201          0.14135
#> 71  profiles5_groups3          5               3   VVI individuals     3     244          0.17159
#> 72  profiles5_groups3          5               3   VVI individuals     4     457          0.32138
#> 73  profiles5_groups3          5               3   VVI individuals     5     397          0.27918
#> 74  profiles5_groups3          5               3   VVI      groups     1      27          0.25472
#> 75  profiles5_groups3          5               3   VVI      groups     2      48          0.45283
#> 76  profiles5_groups3          5               3   VVI      groups     3      31          0.29245
#> 77  profiles6_groups3          6               3   VVI individuals     1     206          0.14487
#> 78  profiles6_groups3          6               3   VVI individuals     2     446          0.31364
#> 79  profiles6_groups3          6               3   VVI individuals     3     322          0.22644
#> 80  profiles6_groups3          6               3   VVI individuals     4     128          0.09001
#> 81  profiles6_groups3          6               3   VVI individuals     5     168          0.11814
#> 82  profiles6_groups3          6               3   VVI individuals     6     152          0.10689
#> 83  profiles6_groups3          6               3   VVI      groups     1      47          0.44340
#> 84  profiles6_groups3          6               3   VVI      groups     2      32          0.30189
#> 85  profiles6_groups3          6               3   VVI      groups     3      27          0.25472
#> 86  profiles2_groups4          2               4   VVI individuals     1     832          0.58509
#> 87  profiles2_groups4          2               4   VVI individuals     2     590          0.41491
#> 88  profiles2_groups4          2               4   VVI      groups     1      25          0.23585
#> 89  profiles2_groups4          2               4   VVI      groups     2      12          0.11321
#> 90  profiles2_groups4          2               4   VVI      groups     3      22          0.20755
#> 91  profiles2_groups4          2               4   VVI      groups     4      47          0.44340
#> 92  profiles3_groups4          3               4   VVI individuals     1     477          0.33544
#> 93  profiles3_groups4          3               4   VVI individuals     2     576          0.40506
#> 94  profiles3_groups4          3               4   VVI individuals     3     369          0.25949
#> 95  profiles3_groups4          3               4   VVI      groups     1      46          0.43396
#> 96  profiles3_groups4          3               4   VVI      groups     2      15          0.14151
#> 97  profiles3_groups4          3               4   VVI      groups     3      20          0.18868
#> 98  profiles3_groups4          3               4   VVI      groups     4      25          0.23585
#> 99  profiles4_groups4          4               4   VVI individuals     1     191          0.13432
#> 100 profiles4_groups4          4               4   VVI individuals     2     321          0.22574
#> 101 profiles4_groups4          4               4   VVI individuals     3     510          0.35865
#> 102 profiles4_groups4          4               4   VVI individuals     4     400          0.28129
#> 103 profiles4_groups4          4               4   VVI      groups     1      46          0.43396
#> 104 profiles4_groups4          4               4   VVI      groups     2      24          0.22642
#> 105 profiles4_groups4          4               4   VVI      groups     3      11          0.10377
#> 106 profiles4_groups4          4               4   VVI      groups     4      25          0.23585
#> 107 profiles5_groups4          5               4   VVI individuals     1     189          0.13291
#> 108 profiles5_groups4          5               4   VVI individuals     2      38          0.02672
#> 109 profiles5_groups4          5               4   VVI individuals     3     302          0.21238
#> 110 profiles5_groups4          5               4   VVI individuals     4     491          0.34529
#> 111 profiles5_groups4          5               4   VVI individuals     5     402          0.28270
#> 112 profiles5_groups4          5               4   VVI      groups     1      23          0.21698
#> 113 profiles5_groups4          5               4   VVI      groups     2      25          0.23585
#> 114 profiles5_groups4          5               4   VVI      groups     3      11          0.10377
#> 115 profiles5_groups4          5               4   VVI      groups     4      47          0.44340
#> 116 profiles6_groups4          6               4   VVI individuals     1     491          0.34529
#> 117 profiles6_groups4          6               4   VVI individuals     2     107          0.07525
#> 118 profiles6_groups4          6               4   VVI individuals     3     341          0.23980
#> 119 profiles6_groups4          6               4   VVI individuals     4     149          0.10478
#> 120 profiles6_groups4          6               4   VVI individuals     5     174          0.12236
#> 121 profiles6_groups4          6               4   VVI individuals     6     160          0.11252
#> 122 profiles6_groups4          6               4   VVI      groups     1      26          0.24528
#> 123 profiles6_groups4          6               4   VVI      groups     2      14          0.13208
#> 124 profiles6_groups4          6               4   VVI      groups     3      23          0.21698
#> 125 profiles6_groups4          6               4   VVI      groups     4      43          0.40566
#>     estimated_n estimated_proportion average_posterior odds_correct_classification relative_entropy
#> 1        832.92              0.58574             0.982                       37.75            0.922
#> 2        589.08              0.41426             0.974                       53.24            0.922
#> 3        452.67              0.31834             0.832                       10.61            0.772
#> 4        551.11              0.38756             0.964                       42.39            0.772
#> 5        418.21              0.29410             0.863                       15.17            0.772
#> 6        478.06              0.33619             0.842                       10.51            0.735
#> 7        172.20              0.12110             0.803                       29.63            0.735
#> 8        405.28              0.28501             0.849                       14.14            0.735
#> 9        366.46              0.25771             0.846                       15.86            0.735
#> 10       359.22              0.25262             0.847                       16.42            0.773
#> 11       405.22              0.28496             0.849                       14.13            0.773
#> 12       480.38              0.33782             0.842                       10.42            0.773
#> 13         4.89              0.00344             0.962                     7368.14            0.773
#> 14       172.30              0.12116             0.804                       29.66            0.773
#> 15       340.87              0.23971             0.735                        8.78            0.669
#> 16       165.78              0.11658             0.784                       27.56            0.669
#> 17       348.85              0.24532             0.780                       10.89            0.669
#> 18        95.07              0.06686             0.636                       24.43            0.669
#> 19       321.57              0.22614             0.812                       14.74            0.669
#> 20       149.86              0.10538             0.627                       14.29            0.669
#> 21       834.10              0.58657             0.985                       46.08            0.938
#> 22       587.90              0.41343             0.982                       76.35            0.938
#> 23        71.52              0.67470             0.992                       59.89            0.940
#> 24        34.48              0.32530             0.969                       64.80            0.940
#> 25       469.56              0.33021             0.849                       11.40            0.790
#> 26       569.85              0.40074             0.979                       68.94            0.790
#> 27       382.59              0.26905             0.850                       15.45            0.790
#> 28        71.25              0.67213             0.992                       62.14            0.943
#> 29        34.75              0.32787             0.977                       87.77            0.943
#> 30       480.19              0.33769             0.845                       10.71            0.739
#> 31       178.19              0.12531             0.791                       26.49            0.739
#> 32       407.65              0.28667             0.858                       15.08            0.739
#> 33       355.96              0.25033             0.846                       16.41            0.739
#> 34        71.16              0.67132             0.990                       48.96            0.937
#> 35        34.84              0.32868             0.975                       80.76            0.937
#> 36       464.18              0.32643             0.846                       11.33            0.775
#> 37       408.11              0.28700             0.860                       15.27            0.775
#> 38        36.89              0.02594             0.838                      194.58            0.775
#> 39       337.67              0.23746             0.848                       17.88            0.775
#> 40       175.15              0.12317             0.803                       28.94            0.775
#> 41        71.40              0.67360             0.990                       48.07            0.965
#> 42        34.60              0.32640             0.996                      581.20            0.965
#> 43       444.31              0.31245             0.835                       11.14            0.732
#> 44       173.96              0.12233             0.781                       25.62            0.732
#> 45       294.72              0.20726             0.728                       10.24            0.732
#> 46       118.27              0.08317             0.699                       25.64            0.732
#> 47       349.82              0.24600             0.845                       16.71            0.732
#> 48        40.93              0.02878             0.789                      126.37            0.732
#> 49        69.10              0.65185             0.987                       39.18            0.946
#> 50        36.90              0.34815             0.972                       65.76            0.946
#> 51       830.91              0.58433             0.986                       49.28            0.940
#> 52       591.09              0.41567             0.982                       75.73            0.940
#> 53        26.66              0.25147             0.905                       28.33            0.816
#> 54        30.78              0.29035             0.962                       62.25            0.816
#> 55        48.57              0.45818             0.900                       10.69            0.816
#> 56       368.95              0.25946             0.845                       15.58            0.792
#> 57       574.96              0.40433             0.978                       65.03            0.792
#> 58       478.08              0.33621             0.853                       11.45            0.792
#> 59        27.65              0.26086             0.898                       24.87            0.797
#> 60        31.03              0.29271             0.977                      104.20            0.797
#> 61        47.32              0.44643             0.864                        7.90            0.797
#> 62       497.03              0.34953             0.854                       10.86            0.739
#> 63       391.99              0.27566             0.847                       14.60            0.739
#> 64       202.20              0.14220             0.801                       24.24            0.739
#> 65       330.78              0.23261             0.842                       17.57            0.739
#> 66        27.40              0.25844             0.874                       19.96            0.823
#> 67        48.01              0.45292             0.915                       13.01            0.823
#> 68        30.60              0.28864             0.981                      127.50            0.823
#> 69       141.26              0.09934             0.628                       15.34            0.714
#> 70       211.64              0.14884             0.804                       23.52            0.714
#> 71       248.08              0.17446             0.804                       19.44            0.714
#> 72       435.18              0.30603             0.817                       10.09            0.714
#> 73       385.84              0.27134             0.843                       14.42            0.714
#> 74        28.08              0.26491             0.917                       30.75            0.854
#> 75        46.73              0.44081             0.918                       14.29            0.854
#> 76        31.19              0.29428             0.987                      180.98            0.854
#> 77       206.70              0.14536             0.777                       20.53            0.674
#> 78       414.48              0.29148             0.795                        9.42            0.674
#> 79       297.05              0.20890             0.733                       10.42            0.674
#> 80       150.21              0.10563             0.711                       20.85            0.674
#> 81       191.71              0.13482             0.655                       12.18            0.674
#> 82       161.85              0.11382             0.788                       28.90            0.674
#> 83        46.26              0.43639             0.925                       15.84            0.853
#> 84        32.68              0.30828             0.984                      137.98            0.853
#> 85        27.07              0.25534             0.915                       31.49            0.853
#> 86       829.91              0.58362             0.986                       50.18            0.942
#> 87       592.09              0.41638             0.984                       85.01            0.942
#> 88        24.59              0.23201             0.868                       21.73            0.798
#> 89        11.33              0.10689             0.851                       47.77            0.798
#> 90        22.34              0.21071             0.877                       26.81            0.798
#> 91        47.74              0.45039             0.912                       12.62            0.798
#> 92       476.54              0.33512             0.853                       11.52            0.793
#> 93       576.30              0.40527             0.980                       70.36            0.793
#> 94       369.16              0.25960             0.843                       15.32            0.793
#> 95        47.04              0.44379             0.902                       11.51            0.772
#> 96        13.62              0.12851             0.841                       35.95            0.772
#> 97        20.96              0.19774             0.868                       26.74            0.772
#> 98        24.37              0.22995             0.843                       17.98            0.772
#> 99       201.57              0.14175             0.801                       24.37            0.740
#> 100      332.98              0.23416             0.843                       17.58            0.740
#> 101      494.01              0.34741             0.850                       10.64            0.740
#> 102      393.44              0.27668             0.853                       15.17            0.740
#> 103       46.89              0.44232             0.903                       11.76            0.803
#> 104       22.90              0.21605             0.867                       23.74            0.803
#> 105       10.91              0.10288             0.898                       76.55            0.803
#> 106       25.31              0.23874             0.873                       21.82            0.803
#> 107      201.07              0.14140             0.803                       24.79            0.751
#> 108       61.45              0.04322             0.786                       81.48            0.751
#> 109      305.33              0.21472             0.831                       17.96            0.751
#> 110      461.34              0.32443             0.826                        9.90            0.751
#> 111      392.80              0.27623             0.853                       15.15            0.751
#> 112       23.02              0.21714             0.926                       44.87            0.830
#> 113       25.69              0.24237             0.882                       23.26            0.830
#> 114       11.14              0.10505             0.899                       75.62            0.830
#> 115       46.16              0.43544             0.909                       12.98            0.830
#> 116      466.75              0.32823             0.830                        9.98            0.725
#> 117      114.80              0.08073             0.847                       63.10            0.725
#> 118      320.60              0.22546             0.806                       14.24            0.725
#> 119      167.55              0.11783             0.765                       24.39            0.725
#> 120      192.57              0.13542             0.797                       25.09            0.725
#> 121      159.73              0.11233             0.654                       14.97            0.725
#> 122       25.05              0.23630             0.898                       28.53            0.851
#> 123       13.24              0.12488             0.907                       68.41            0.851
#> 124       23.33              0.22010             0.919                       40.18            0.851
#> 125       44.38              0.41872             0.940                       21.71            0.851
```

`plot(x, what = "diagnostics")` draws relative entropy and the lowest
AvePP across the candidates, with the cutoff dashed; `what = "entropy"`
and `what = "avepp"` draw one of the two. For one fitted candidate,
[`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
draws the certainty of every unit in its class, the uncertainty of every
unit, the class sizes and the average posterior matrix, for the profiles
and, when there are two or more group classes, for the group classes
(`what = "group_posteriors"` and the like). The two levels often differ
sharply: a group is classified from all of its observations together,
and an observation from its own responses alone.

``` r

plot(models, what = "diagnostics")
```

![](model-selection_files/figure-html/diagnostics-1.png)

Entropy and AvePP can disagree, and the difference is informative. AvePP
reads only the posterior of each unit’s assigned class; entropy reads
the whole row of posteriors. A unit at 0.70 for its class and 0.10 for
each of three others, and a unit at 0.70 and 0.30, give the same AvePP,
but the first, spread over every class, lowers entropy much more.
Entropy also falls faster for a few very uncertain units than it rises
for many certain ones. A level with most units near certainty and a tail
of uncertain ones therefore shows a moderate AvePP beside a low entropy.

The two levels of a two-level model are classified from different
amounts of information. A group class is assigned from all of a group’s
observations together, an observation’s profile from its own responses
alone, so group classes are usually much sharper. For one candidate,
[`plot_diagnostics()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
draws both levels:

``` r

plot_diagnostics(get_candidate(models, "profiles4_groups3"),
                 what = "group_posteriors")
```

![](model-selection_files/figure-html/group-certainty-1.png)

When every candidate separates its profiles poorly, so that no number of
classes brings entropy or the smallest AvePP up, the indicators do not
separate single observations well, and the number of classes is chosen
from the criteria and the interpretation, with the separation reported
as a limitation. Modal assignments are then a poor substitute for the
classes in a later analysis, and the bias-corrected three-step estimates
([`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md),
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md))
are used instead.

Classification quality is not a criterion for the number of classes in
itself: a model with fewer classes is usually better separated, so
choosing by entropy alone favours too few classes. It decides between
candidates that the criteria leave close, and it qualifies how the
chosen classes can be used, for instance whether modal assignments can
stand in for the classes in a later analysis.

## Can the estimate be trusted?

### Estimation checks

The likelihood of a mixture model has several local maxima. A solution
reached from several random starts is more credible than one reached
once, and an estimate on the boundary of the parameter space, a variance
at its floor or a class with almost no members, marks a model too
complex for the data. `plot(x, what = "trust")` gives each candidate one
tile, labelled with the number of starts that reached its best log
likelihood, and coloured by the first problem it has:

| Status | Meaning |
|----|----|
| best replicated | converged, inside the parameter space, best likelihood reached by more than one start |
| best found once | as above, but only one start reached the best likelihood |
| on a boundary | an estimate reached the edge of the parameter space |
| not converged | the EM algorithm stopped before meeting its tolerance |
| not fitted (`x`) | not identified, or failed with an error kept in the table |

``` r

plot(models, what = "trust")
```

![](model-selection_files/figure-html/trust-1.png)

A candidate found once is refitted with more starts before it is chosen,
and so is the rest of the grid (see *Random starts* above). If the refit
raises its likelihood, the earlier solution was a local maximum and the
new one replaces it; its entropy and AvePP often rise with it. A higher
AvePP is never a reason to prefer one solution over another of the same
candidate: the solution with the higher likelihood is the estimate. The
columns `converged`, `boundary`, `n_best_replicated`, `warnings` and
`error` of the table of candidates hold the same checks.

### The smallest class

A class holding a very small share of the units is hard to interpret and
estimated from little data. `plot(x, what = "smallest")` draws the
estimated share of the smallest class of every candidate, with a dashed
line at 5%, a common floor rather than a rule.

``` r

plot(models, what = "smallest")
```

![](model-selection_files/figure-html/smallest-1.png)

## All of the plots together

`plot(x, what = "selection")` returns the deltas, the heatmap, the
smallest class, the estimation checks and the classification diagnostics
as one list of plots, drawn in turn. `using` sets the criterion of the
first two and `cutoff` the AvePP line of the last.

``` r

plot(models, what = "selection")
```

![](model-selection_files/figure-html/selection-1.png)![](model-selection_files/figure-html/selection-2.png)![](model-selection_files/figure-html/selection-3.png)![](model-selection_files/figure-html/selection-4.png)![](model-selection_files/figure-html/selection-5.png)

## How the profiles divide

The tree shows how the profiles of each model divide into those of the
next: each profile is linked to the profiles of the model with one more
that receive its observations (Zappia and Oshlack, 2018). A profile that
splits cleanly in two suggests a real distinction; one that draws from
several profiles of the smaller model suggests a reshuffling.

``` r

plot(models, what = "tree")
```

![](model-selection_files/figure-html/tree-1.png)

## Comparing chosen candidates

Once the evidence narrows the choice to a few candidates,
[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
places them side by side, named or by their numbers of classes. Its
table adds the difference of each candidate’s BIC from the lowest among
those compared, `delta_bic`, and its BIC weight, `exp(-delta_bic / 2)`
normalized over the compared candidates: the approximate posterior
probability of each if one of them is true and all were equally likely
beforehand (Wagenmakers and Farrell, 2004). The weights depend on which
candidates are compared.

``` r

comparison <- compare_models(models, "profiles3_groups3", "profiles4_groups3",
                             "profiles4_groups2")
comparison
#> Model comparison
#> 
#> Model              Profiles  Group classes   k    LogLik       BIC   ΔBIC  Weight       ICL  Entropy  Min AvePP  < 0.8  Group entropy  Group min AvePP
#> -----------------  --------  -------------  --  --------  --------  -----  ------  --------  -------  ---------  -----  -------------  ---------------
#> profiles3_groups3         3              3  38  -8348.43  16874.07  21.92   0.000  17622.03     0.79       0.85      0           0.80             0.86
#> profiles4_groups3         4              3  51  -8307.15  16852.14   0.00   1.000  18014.13     0.74       0.80      0           0.82             0.87
#> profiles4_groups2         4              2  47  -8326.19  16871.57  19.43   0.000  18021.23     0.74       0.79      1           0.94             0.98
plot(comparison)
```

![](model-selection_files/figure-html/compare-1.png)

The model behind a row is taken out with
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md),
by name or by numbers of classes.

``` r

chosen <- get_candidate(models, "profiles4_groups3")
chosen
#> Two-level latent profile analysis: 4 profiles, 3 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -8307.154320 | AIC: 16716.309 | BIC (groups): 16852.144
#> Converged: TRUE | iterations: 48 | best start: 3/3
#> 
#>  profile browse lectures forum_read forum_post attendance count proportion
#>        1  0.412    0.232      0.490      0.383      0.309   497      0.350
#>        2 -0.564   -0.561     -0.695     -0.676     -0.656   392      0.276
#>        3 -1.113   -0.685     -1.200     -0.880     -1.417   202      0.142
#>        4  0.728    0.734      0.821      0.763      1.180   331      0.233
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

## Testing one more class

The usual likelihood ratio test is invalid between mixtures with `k` and
`k + 1` classes, because the smaller model sits on the boundary of the
larger one, where the statistic loses its reference distribution.
[`bootstrap_lrt()`](https://pak.dynasite.org/latents/reference/bootstrap_lrt.md)
simulates that distribution (McLachlan, 1987; Nylund, Asparouhov and
Muthén, 2007): data sets are drawn from the smaller fitted model, both
models are refitted to each, and the observed statistic is placed among
the simulated ones. With `B` simulated data sets and `b` statistics at
least as large as the observed one, the p-value is `(b + 1) / (B + 1)`,
which can never be zero and cannot fall below `1 / (B + 1)`. A replicate
that fails a fitting check withholds the p-value rather than drop
silently. The example uses 19 data sets, whose smallest attainable
p-value is 0.05; applied analyses use 99 or more.

``` r

smaller <- get_candidate(models, "profiles3_groups3")
larger <- get_candidate(models, "profiles4_groups3")
test <- bootstrap_lrt(smaller, larger, data = course_engagement,
                      iter = 19, n_starts = 2, seed = 1)
test
#> Parametric bootstrap likelihood-ratio comparison
#> Null: 3 profiles, 3 group classes; alternative: 4 profiles, 3 group classes
#> Observed statistic: 82.548767
#> p-value: 0.05 (Monte Carlo SE 0.0487) from 19 of 19 valid replicates
plot(test)
```

![](model-selection_files/figure-html/blrt-1.png)

## Putting the evidence together

No rule turns this evidence into a number of classes, and latents does
not apply one. A choice is made by reading the evidence in a fixed order
and reporting every step:

0.  **State the preparation**: the number of starts, any centring and
    scaling, and how missing values were handled, with how many rows had
    nothing observed.
1.  **State the convention** for the sample size before reading any
    criterion, and keep it.
2.  **Read the minima** of the criteria under that convention. Note
    where they agree and in which direction they disagree. A minimum on
    the edge of the grid calls for a wider grid.
3.  **Read the deltas and knees** when the criteria keep falling. Report
    the knee with its `distance` and whether it is `stable`.
4.  **Check classification** of the candidates still in contention:
    relative entropy, the smallest AvePP and OCC of each level.
5.  **Check estimation**: replication of the best likelihood, boundary
    estimates and the smallest class. Refit a candidate found once with
    more starts.
6.  **Test** the leading candidate against its neighbour with the
    bootstrap likelihood ratio test when the evidence is close.
7.  **Interpret** the leading candidates. Between candidates that the
    statistics leave close, the one whose classes are distinct and
    meaningful for the question is chosen, and the reason is stated.

The tables behind every step are returned by the verbs above, so the
reported choice can be reproduced from the reported evidence.

## Functions

| Call | Result |
|:---|:---|
| enumerate_classes(…, n_starts = , centering = , scaling = , missing = “fiml”) | The same starts, centring, scaling and missing-data handling for every candidate |
| as.data.frame(models) | One row per candidate: every criterion, classification summaries and checks |
| get_results(models, “best”) | The candidate at the minimum of each criterion and convention |
| plot(models), plot(models, what = “bic_groups”), what = “all” | Criteria across the candidates, lowest ringed |
| get_results(models, “deltas”, using = ) | Decrease in a criterion from one profile fewer and one group class fewer |
| plot(models, what = “deltas”, using = ) | The decreases along profiles and along group classes |
| get_results(models, “knees”, using = ) | Chord-distance knee of each line, with its stability |
| plot(models, what = “heatmap”, using = ) | Distance of every candidate from the lowest criterion |
| get_results(models, “classification”) | Size, AvePP, OCC and entropy of every class of every candidate |
| get_results(models, “candidates”, cutoff = ) | The table of candidates with another AvePP cutoff |
| plot(models, what = “diagnostics”, “entropy” or “avepp”, cutoff = ) | Relative entropy and lowest AvePP across the candidates |
| get_results(fit, “entropy”) | Relative entropy of one fit at each level, over the classified units |
| plot_diagnostics(fit), what = “group_posteriors” and the like | Certainty, uncertainty, sizes and AvePP of one fit, profiles and group classes |
| plot(models, what = “trust”) | Convergence, boundary and replication of every candidate |
| plot(models, what = “smallest”) | Share of the smallest class of every candidate |
| plot(models, what = “selection”) | Deltas, heatmap, smallest class, trust and diagnostics together |
| plot(models, what = “tree”) | How the profiles divide as more are added |
| compare_models(models, “name”, …) | Chosen candidates side by side, with BIC weights |
| get_candidate(models, “name”) | One fitted candidate |
| bootstrap_lrt(smaller, larger, data, iter) | Bootstrap likelihood ratio test of k against k + 1 classes |

## References

Akaike, H. (1974). A new look at the statistical model identification.
*IEEE Transactions on Automatic Control*, 19, 716–723.

Banfield, J. D., & Raftery, A. E. (1993). Model-based Gaussian and
non-Gaussian clustering. *Biometrics*, 49, 803–821.

Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture
model for clustering with the integrated completed likelihood. *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 22, 719–725.

Biernacki, C., & Govaert, G. (1997). Using the classification likelihood
to choose the number of clusters. *Computing Science and Statistics*,
29, 451–457.

Bozdogan, H. (1987). Model selection and Akaike’s information criterion.
*Psychometrika*, 52, 345–370.

Cavanaugh, J. E. (1999). A large-sample model selection criterion based
on Kullback’s symmetric divergence. *Statistics and Probability
Letters*, 42, 333–343.

McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test
statistic for the number of components in a normal mixture. *Applied
Statistics*, 36, 318–324.

Nagin, D. S. (2005). *Group-based modeling of development*. Harvard
University Press.

Nylund, K. L., Asparouhov, T., & Muthén, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, 14,
535–569.

Satopää, V., Albrecht, J., Irwin, D., & Raghavan, B. (2011). Finding a
“kneedle” in a haystack: detecting knee points in system behavior. *31st
International Conference on Distributed Computing Systems Workshops*,
166–171.

Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics*, 6, 461–464.

Sclove, S. L. (1987). Application of model-selection criteria to some
problems in multivariate analysis. *Psychometrika*, 52, 333–343.

Wagenmakers, E.-J., & Farrell, S. (2004). AIC model selection using
Akaike weights. *Psychonomic Bulletin & Review*, 11, 192–196.

Zappia, L., & Oshlack, A. (2018). Clustering trees: a visualization for
evaluating clusterings at multiple resolutions. *GigaScience*, 7(7),
giy083.
