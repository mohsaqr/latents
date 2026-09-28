# Mixture regression

A single regression assumes that one line describes everyone. Mixture
regression drops that assumption: the data are taken to come from a
small number of latent classes, each with its own regression, and the
model estimates the regressions, how common each class is, and how
probably each observation belongs to each class — all at once. The
method goes by several names: clusterwise regression, latent class
regression, regression mixture (DeSarbo and Cron 1988; Wedel and DeSarbo
1995).

[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
fits it for continuous (`"gaussian"`), binary or binomial (`"binomial"`)
and count (`"poisson"`) outcomes, at one level or two:

| Question | Call |
|----|----|
| Do observations follow different regressions? | `mixture_regression(y ~ x, data, n_classes)` |
| Do *groups* (persons, schools) follow different regressions? | `mixture_regression(..., id = "person", class_level = "group")` |
| Do groups differ in how often their observations follow each regression? | `mixture_regression(..., id = "person", n_group_classes = 2)` |
| What predicts which class an observation (or group) is in? | `membership = ~ z`, `group_membership = ~ v` |
| Is some effect the same in every class? | `common = ~ x2` |
| How many classes? | `enumerate_regressions(..., n_classes = 1:4)` |

## The data

`study_hours` is simulated so that the truth is known. 150 students
report six weeks each. In each week a student studies with a *deep* or a
*surface* strategy, and the strategy decides how hours turn into a quiz
score: deep weeks gain about 4.5 points per hour from a base of 35,
surface weeks under a point per hour from a base of 55. Steady students
study deeply in most weeks, erratic students in few. The generating
strategy and student type are kept in the data so recovery can be
checked; the models never see them.

``` r

library(latents)
head(study_hours)
#>   student week hours sleep motivation score passed questions strategy student_type
#> 1       1    1  10.2   5.7      -0.08  67.8      0         2  surface      erratic
#> 2       1    2   2.6   6.4      -0.08  59.5      0         2  surface      erratic
#> 3       1    3  10.1   6.4      -0.08  55.9      1         2  surface      erratic
#> 4       1    4   8.7   6.7      -0.08  52.5      0         5  surface      erratic
#> 5       1    5  11.3   9.5      -0.08  89.4      0         2     deep      erratic
#> 6       1    6   5.7   6.7      -0.08  58.0      0         3  surface      erratic
```

## One class per observation

``` r

fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 5, seed = 1)
fit
#> Mixture of gaussian regressions: 2 classes (one class per row)
#> score ~ hours
#> 900 rows | log likelihood -3144.0074 | BIC 6335.63 | entropy 0.443
#> Converged: TRUE | iterations: 42 | best likelihood reached by 5 of 5 completed starts (6 run)
#> 
#>    class        term estimate std_error   p_value
#>  class_1 (Intercept)  34.3545    0.7193 0.000e+00
#>  class_1       hours   4.5659    0.1047 0.000e+00
#>  class_2 (Intercept)  54.8951    1.0927 0.000e+00
#>  class_2       hours   0.8084    0.1696 1.876e-06
#> 
#> Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".
```

The two classes are the two strategies: a steep regression from a low
base and a flat one from a high base. Every table comes from
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md);
the class table gives each class’s share, its residual standard
deviation and how certainly its members are assigned.

``` r

get_results(fit, "classes")
#>     class share count n_assigned mean_posterior sigma sigma_std_error prior prior_std_error
#> 1 class_1 0.564   508        573          0.788  5.24           0.278 0.564           0.035
#> 2 class_2 0.436   392        327          0.829  6.89           0.346 0.436           0.035
```

`get_results(fit, "assignments")` gives one row per observation. Because
these data carry the generating strategy, the recovery table
cross-tabulates the assigned classes against it:

``` r

get_results(fit, "recovery", data = study_hours, truth = "strategy")
#>   assigned strategy   n share
#> 1  class_1     deep 448 0.782
#> 2  class_2     deep  49 0.150
#> 3  class_1  surface 125 0.218
#> 4  class_2  surface 278 0.850
```

``` r

plot(fit)
```

![](mixture-regression_files/figure-html/plot-fitted-1.png)

## How many classes?

[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md)
fits a range of class counts and returns their fit statistics in one
table, marking the minimum-BIC model.

``` r

classes <- enumerate_regressions(score ~ hours, data = study_hours,
                                 n_classes = 1:4, n_starts = 3, seed = 1)
classes
#> Mixture regression class enumeration
#> 
#>  n_classes n_group_classes log_likelihood n_parameters  bic  icl entropy smallest_share converged
#>          1               1          -3245            3 6510 6510      NA        1.00000      TRUE
#>          2               1          -3144            7 6336 7030  0.4432        0.43526      TRUE
#>          3               1          -3133           11 6341 7272  0.5292        0.07523      TRUE
#>          4               1          -3126           15 6354 7460  0.5571        0.07192      TRUE
#>  best_bic
#>     FALSE
#>      TRUE
#>     FALSE
#>     FALSE
#> 
#> One model: get_results(x, "model", n_classes = ).
```

With `bootstrap = B`, each row also gets a parametric bootstrap
likelihood-ratio test of one class fewer (McLachlan 1987): outcomes are
simulated from the smaller model, both models are refitted, and the
p-value is the share of simulated statistics at least as large as the
observed one. It is slow — `2 B` refits per row — so it is not run here.

## What predicts the class?

`membership` adds a multinomial logit for class membership. Here sleep
before the quiz makes a deep week more likely, so it should lower the
odds of the surface class:

``` r

with_sleep <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                                 membership = ~ sleep, n_starts = 5, seed = 1)
get_results(with_sleep, "membership")
#>   model   class        term reference estimate std_error statistic p_value conf_low conf_high
#> 1 class class_2 (Intercept)   class_1    1.555    0.6817      2.28 0.02255    0.219    2.8911
#> 2 class class_2       sleep   class_1   -0.261    0.0962     -2.72 0.00661   -0.450   -0.0727
#>   odds_ratio
#> 1       4.74
#> 2       0.77
```

The coefficients are estimated jointly with the regressions — a one-step
model, not a regression on assigned classes — so classification error
does not attenuate them.

## Groups that share a class

When every row of a group belongs to the same class, pass `id` and
`class_level = "group"`. A group’s evidence for a class is the product
of its rows’ likelihoods, so classes of *students* are separated far
more sharply than classes of weeks could be. The binary `passed` outcome
in these data depends on the student’s type, not on the week:

``` r

students <- mixture_regression(passed ~ hours, data = study_hours, n_classes = 2,
                               family = "binomial", id = "student",
                               class_level = "group", membership = ~ motivation,
                               n_starts = 5, seed = 1)
students
#> Mixture of binomial regressions: 2 classes (one class per `student`)
#> passed ~ hours
#> 900 rows in 150 groups | log likelihood -540.8146 | BIC 1111.69 | entropy 0.668
#> Converged: TRUE | iterations: 16 | best likelihood reached by 6 of 6 completed starts (6 run)
#> 
#>    class        term estimate std_error   p_value
#>  class_1 (Intercept) -2.46252   0.32414 3.031e-14
#>  class_1       hours  0.51369   0.06253 2.118e-16
#>  class_2 (Intercept)  0.07049   0.26795 7.925e-01
#>  class_2       hours -0.04109   0.04113 3.178e-01
#> 
#> Every table: get_results(x, what = ), e.g. "classes", "membership", "fit", "assignments".
get_results(students, "membership")
#>   model   class        term reference estimate std_error statistic  p_value conf_low conf_high
#> 1 class class_2 (Intercept)   class_1    -0.15     0.374    -0.402 0.687587   -0.883     0.583
#> 2 class class_2  motivation   class_1    -1.82     0.491    -3.703 0.000213   -2.781    -0.856
#>   odds_ratio
#> 1      0.860
#> 2      0.162
```

A mixture of logistic regressions with a single binary outcome per class
assignment is not identified — any mixture of Bernoulli distributions is
a Bernoulli distribution (Follmann and Lambert 1991) — so
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
refuses `family = "binomial"` with one trial per row unless the class is
shared by a group:

``` r

mixture_regression(passed ~ hours, data = study_hours, n_classes = 2,
                   family = "binomial")
#> Error:
#> ! A mixture of logistic regressions with one binary outcome per class assignment is not identified: a mixture of Bernoulli distributions is itself a Bernoulli distribution (Follmann and Lambert 1991). Give each class assignment several outcomes with `id` and `class_level = "group"`, or supply binomial counts as `cbind(successes, failures)`.
```

## Two levels: groups that differ in their mix

The deep/surface strategy changes from week to week, but steady and
erratic students use the strategies in different proportions. That is a
two-level mixture (Vermunt 2003): rows keep their own regression class,
and a latent *group class* shifts the class probabilities of all rows in
its groups.

``` r

two_level <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                                id = "student", n_group_classes = 2,
                                membership = ~ sleep, group_membership = ~ motivation,
                                n_starts = 5, seed = 1)
summary(two_level)
#> Model fit
#>    family   nesting n_classes n_group_classes n_observations n_groups log_likelihood n_parameters
#>  gaussian two-level         2               2            900      150          -3092           11
#>   aic  bic bic_rows sabic  icl entropy group_entropy smallest_share converged iterations n_starts
#>  6206 6239     6258  6204 6785  0.5621        0.6646         0.4207      TRUE         53        6
#>  n_best_replicated vcov_type
#>                  3  observed
#> 
#> Regression coefficients
#>    class        term estimate std_error statistic   p_value conf_low conf_high p_adjusted
#>  class_1 (Intercept)  34.6314   0.61084    56.694 0.000e+00  33.4342    35.829         NA
#>  class_1       hours   4.5230   0.08418    53.728 0.000e+00   4.3580     4.688  0.000e+00
#>  class_2 (Intercept)  55.3412   0.90779    60.963 0.000e+00  53.5619    57.120         NA
#>  class_2       hours   0.7182   0.13666     5.255 1.478e-07   0.4503     0.986  1.478e-07
#> 
#> Classes
#>    class  share count n_assigned mean_posterior sigma sigma_std_error
#>  class_1 0.5793 521.4        532         0.8732 5.281          0.2232
#>  class_2 0.4207 378.6        368         0.8455 6.789          0.3134
#> 
#> Class membership (multinomial logit)
#>        model         class                      term     reference estimate std_error statistic
#>  group_class group_class_2               (Intercept) group_class_1  -0.1257    0.4572    -0.275
#>  group_class group_class_2                motivation group_class_1  -1.8663    0.4760    -3.921
#>        class       class_2 (Intercept):group_class_1       class_1   1.3446    0.9490     1.417
#>        class       class_2 (Intercept):group_class_2       class_1   4.2492    1.0602     4.008
#>        class       class_2                     sleep       class_1  -0.4688    0.1392    -3.368
#>    p_value conf_low conf_high odds_ratio
#>  7.833e-01  -1.0219    0.7704     0.8818
#>  8.834e-05  -2.7992   -0.9333     0.1547
#>  1.565e-01  -0.5155    3.2046     3.8365
#>  6.122e-05   2.1713    6.3271    70.0498
#>  7.566e-04  -0.7416   -0.1960     0.6257
#> 
#> Group classes
#>    group_class   class probability group_share n_groups_assigned
#>  group_class_1 class_1      0.8609      0.5134                76
#>  group_class_1 class_2      0.1391      0.5134                76
#>  group_class_2 class_1      0.2821      0.4866                74
#>  group_class_2 class_2      0.7179      0.4866                74
```

The `group_classes` table is the heart of it: within each group class,
the probability of each regression class. The group-class membership
model shows motivation separating the two kinds of student, and the
class model shows the week-level effect of sleep.

``` r

get_results(two_level, "recovery", data = study_hours,
            truth = "student_type", by = "group_class")
#>        assigned student_type  n  share
#> 1 group_class_1       steady 69 0.9079
#> 2 group_class_2       steady  8 0.1081
#> 3 group_class_1      erratic  7 0.0921
#> 4 group_class_2      erratic 66 0.8919
```

## Counts, shared effects and equal variances

The Poisson family takes counts, and
[`offset()`](https://rdrr.io/r/stats/offset.html) terms for exposure.
`common` names predictors whose coefficient is the same in every class;
for the Gaussian family `variance = "equal"` shares one residual
standard deviation.

``` r

asked <- mixture_regression(questions ~ hours + sleep, data = study_hours, n_classes = 2,
                            family = "poisson", common = ~ sleep, n_starts = 5, seed = 1)
get_results(asked, "coefficients")
#>     class        term estimate std_error statistic  p_value conf_low conf_high p_adjusted
#> 1 class_1 (Intercept)   0.5578    0.1606     3.474 5.13e-04   0.2431    0.8726         NA
#> 2 class_1       hours   0.1038    0.0107     9.666 4.19e-22   0.0828    0.1248   1.26e-21
#> 3 class_2 (Intercept)   1.5081    0.2490     6.057 1.39e-09   1.0201    1.9962         NA
#> 4 class_2       hours  -0.1099    0.0446    -2.462 1.38e-02  -0.1973   -0.0224   2.07e-02
#> 5  common       sleep  -0.0194    0.0213    -0.914 3.61e-01  -0.0611    0.0222   3.61e-01
#>   exp_estimate
#> 1        1.747
#> 2        1.109
#> 3        4.518
#> 4        0.896
#> 5        0.981
```

`exp_estimate` is the rate ratio of each coefficient.

## Uncertainty

Standard errors come from the observed information by default. Scores
are analytic (the posterior expectation of the complete-data score), and
the information is their numerical derivative. `vcov_type = "robust"`
gives the sandwich estimator, clustered on the independent unit — the
row, or the group when `id` defines one, including a single-level fit
given `id`.

``` r

get_results(fit, "coefficients", vcov_type = "robust")
#>     class        term estimate std_error statistic  p_value conf_low conf_high p_adjusted
#> 1 class_1 (Intercept)   34.355     0.810     42.40 0.00e+00   32.767     35.94         NA
#> 2 class_1       hours    4.566     0.114     39.96 0.00e+00    4.342      4.79   0.00e+00
#> 3 class_2 (Intercept)   54.895     1.117     49.13 0.00e+00   52.705     57.08         NA
#> 4 class_2       hours    0.808     0.176      4.58 4.56e-06    0.463      1.15   4.56e-06
plot(fit, what = "coefficients")
```

![](mixture-regression_files/figure-html/robust-1.png)

`p_adjusted` applies the Benjamini-Hochberg correction across the
slopes.

## Prediction and simulation

[`predict()`](https://rdrr.io/r/stats/predict.html) gives the
class-averaged mean (`"response"`), every class’s mean
(`"class_response"`), or — when the outcome is present — the posterior
class probabilities of new rows (`"posterior"`).
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new outcomes
from the fitted model, which is what the bootstrap test uses.

``` r

new_weeks <- data.frame(hours = c(1, 6, 11))
predict(fit, new_weeks, type = "class_response")
#>   row   class prior fitted
#> 1   1 class_1 0.564   38.9
#> 2   2 class_1 0.564   61.7
#> 3   3 class_1 0.564   84.6
#> 4   1 class_2 0.436   55.7
#> 5   2 class_2 0.436   59.7
#> 6   3 class_2 0.436   63.8
```

## Checking a solution

Mixture likelihoods have local maxima.
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
starts EM from the residuals of a pooled regression and from `n_starts`
random partitions, and the `starts` table shows how many starts reached
the best likelihood. A solution reached by one start only deserves more
starts.

``` r

get_results(fit, "starts")
#>   start      kind log_likelihood converged iterations     stage degenerate selected
#> 1     1 residuals          -3144      TRUE         49 completed      FALSE    FALSE
#> 2     2    random          -3144      TRUE         47 completed      FALSE    FALSE
#> 3     3    random          -3144     FALSE         50  screened      FALSE    FALSE
#> 4     4    random          -3144      TRUE         42 completed      FALSE     TRUE
#> 5     5    random          -3144      TRUE         46 completed      FALSE    FALSE
#> 6     6    random          -3144      TRUE         43 completed      FALSE    FALSE
plot(fit, what = "posteriors")
```

![](mixture-regression_files/figure-html/starts-1.png)

## How this compares with other software

The single-level and group-level models match `flexmix` (Leisch 2004):
the same likelihood at the same parameters, and the same binomial and
Poisson estimates. `flexmix` divides the Gaussian residual sum of
squares by `n - p`, so its Gaussian estimates are not quite the
maximum-likelihood ones and
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
reaches a slightly higher likelihood. `flexmix` has no two-level model,
no analytic standard errors (it refits with
[`optim()`](https://rdrr.io/r/stats/optim.html)), and no bootstrap
likelihood-ratio test.

## References

DeSarbo, W. S., & Cron, W. L. (1988). A maximum likelihood methodology
for clusterwise linear regression. *Journal of Classification*, 5,
249–282.

Follmann, D. A., & Lambert, D. (1991). Identifiability of finite
mixtures of logistic regression models. *Journal of Statistical Planning
and Inference*, 27, 375–381.

Leisch, F. (2004). FlexMix: A general framework for finite mixture
models and latent class regression in R. *Journal of Statistical
Software*, 11(8).

McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test
statistic for the number of components in a normal mixture. *Applied
Statistics*, 36, 318–324.

Vermunt, J. K. (2003). Multilevel latent class models. *Sociological
Methodology*, 33, 213–239.

Wedel, M., & DeSarbo, W. S. (1995). A mixture likelihood approach for
generalized linear models. *Journal of Classification*, 12, 21–55.
