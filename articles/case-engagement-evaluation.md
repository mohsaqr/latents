# Case study: Evaluating multilevel latent profile models

Multilevel latent profile analysis (MLPA) estimates latent profiles at
the observation level while simultaneously modeling how the
probabilities of those profiles vary across higher-level latent classes.
In this framework, the observation-level profiles are defined by their
characteristic indicator means, whereas the higher-level classes are
distinguished by the probabilities with which those profiles occur.

Model evaluation therefore involves several related questions. First, it
is important to determine whether the estimated profiles provide
meaningful separation among observations. Second, the measurement
assumptions should be examined to assess whether the model adequately
represents the relationships among the observed indicators. Third,
alternative model specifications can be compared to determine whether a
different number of observation-level profiles or higher-level classes
provides a more appropriate representation of the data.

The three-profile, three-class model introduced in [the introductory
vignette](https://pak.dynasite.org/latents/articles/case-engagement-profiles.md)
serves as the starting point for these evaluation procedures. The
following sections use this model as a reference while examining profile
separation, measurement assumptions, and alternative model
specifications.

``` r

library(latents)
```

## Data and starting model

The `course_engagement` dataset contains repeated course enrolments
nested within students. The bundled data are synthetic, and their
relationship to the longitudinal engagement study is described in the
[introductory
vignette](https://pak.dynasite.org/latents/articles/case-engagement-profiles.md).
Engagement is measured using five course-level indicators: course-page
views (`browse`), lecture-video views (`lectures`), forum posts read
(`forum_read`), forum posts written (`forum_post`), and days active
(`attendance`). Each indicator is simulated as a continuous value
representing a log-count scale and standardized within course; no raw
counts or [`log1p()`](https://rdrr.io/r/base/Log.html) transformation
are generated. The vector `vars` holds their names.

``` r

vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
```

The variable `student` identifies the higher-level unit, while
`sequence` records the order of each student’s course enrolments. In the
present model, enrolments are treated as observations nested within
students. Courses themselves are not modeled as an additional grouping
level. The multilevel structure therefore distinguishes between
variation across enrolments and variation across students.

The starting model contains three Gaussian enrolment-level profiles and
three student-level latent classes. A central quantity in the model is
the conditional probability

``` math
\pi_{k\mid m}=P(C_{it}=k\mid G_i=m),
```

where $`C_{it}`$ denotes the latent profile of student $`i`$ at occasion
$`t`$, and $`G_i`$ denotes that student’s higher-level latent class.
Each probability lies between zero and one, and the profile
probabilities sum to one within each student class. These probabilities
describe how frequently the different engagement profiles are expected
to occur among students belonging to each latent class. The
profile-specific means and variances are shared across student classes,
so the classes differ in their profile probabilities rather than in the
definitions of the profiles themselves.

The model uses a diagonal covariance structure, which assumes that the
engagement indicators are conditionally independent within each latent
profile. It also assumes that enrolments are conditionally independent
given the student’s latent class and that students are independent of
one another. In practice, these assumptions may be violated. For
example, indicators may remain associated because several measures
reflect related forms of online activity, engagement may persist across
adjacent courses taken by the same student, or students may share
unmodeled course-level influences.

Like other finite-mixture models, MLPA can also be sensitive to local
likelihood maxima. Using multiple starting values reduces the risk that
estimation will settle on a suboptimal local solution, although it does
not guarantee that the global maximum has been identified. Model
evaluation should therefore consider both the substantive interpretation
of the solution and evidence concerning the stability of the likelihood
across different starting values.

The
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
function estimates the model by maximum likelihood. A stricter
convergence tolerance (`tol`) is used here so that the optimization
proceeds closer to a stationary solution before standard errors are
calculated. This provides a more stringent convergence criterion for the
starting model used in the subsequent evaluation steps.

``` r

fit <- multilpa(course_engagement, vars, id = "student", time = "sequence",
                n_profiles = 3, n_group_classes = 3, n_starts = 10,
                max_iter = 2000, tol = 1e-10, seed = 1)
fit
#> Two-level latent profile analysis: 3 profiles, 3 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -8348.428704 | AIC: 16772.857 | BIC (groups): 16874.068
#> Converged: TRUE | iterations: 45 | best log likelihood reached by 10 of 10 starts (10 converged)
#> Smallest standardized Gaussian eigenvalue: 0.2003 (profiles, profile_1); Gaussian mass: 478.1; weight effective n: 603.8.
#> 
#> Class sizes   profile_1 33.6% (478)   profile_2 40.4% (575)   profile_3 25.9% (369)
#> Group classes group_class_1 29.3% (31)   group_class_2 26.1% (28)   group_class_3 44.6% (47)
#> 
#> Means
#>   indicator profile_1 profile_2 profile_3
#>      browse     0.387    -0.777     0.709
#>    lectures     0.186    -0.612     0.713
#>  forum_read     0.446    -0.891     0.811
#>  forum_post     0.338    -0.762     0.750
#>  attendance     0.256    -0.940     1.133
#> 
#> Variances and standard errors: get_results(x, "profiles").
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

The fitted model summarizes 1,422 course enrolments contributed by 106
students. The three estimated profiles show a clear ordering across the
engagement indicators. Profile 3 has the highest mean on every
indicator, Profile 1 has intermediate means, and Profile 2 has the
lowest means. For example, the estimated attendance means are 1.133,
0.256, and -0.940 standard deviations relative to the corresponding
course mean for Profiles 3, 1, and 2, respectively.

These differences in indicator means define the substantive
interpretation of the fitted profiles. In this example, the profiles can
therefore be understood as representing progressively lower levels of
engagement across the observed activity measures. However, the presence
of distinct estimated means does not by itself establish that the
profiles are well supported. The diagnostic checks that follow examine
the extent to which the data support this separation and whether the
fitted solution is numerically stable.

The command `get_results(fit, "starts")` summarizes the outcome of each
model initialization, including convergence information and the
maximized log likelihood. When multiple converged starts arrive at the
same or very similar log-likelihood value, this provides evidence that
the selected solution is numerically stable. By contrast, substantially
different likelihood values across converged starts may indicate the
presence of competing local maxima. Similarly, variance estimates that
reach active parameter bounds can signal estimation difficulties or an
inadequately specified solution. Either pattern warrants further
investigation before the fitted profiles are interpreted in detail.

``` r

get_results(fit, "starts")
#>    start log_likelihood converged iterations error boundary
#> 1      1       -8348.43      TRUE         45  <NA>    FALSE
#> 2      2       -8348.43      TRUE         42  <NA>    FALSE
#> 3      3       -8348.43      TRUE         45  <NA>    FALSE
#> 4      4       -8348.43      TRUE         39  <NA>    FALSE
#> 5      5       -8348.43      TRUE         42  <NA>    FALSE
#> 6      6       -8348.43      TRUE         36  <NA>    FALSE
#> 7      7       -8348.43      TRUE         36  <NA>    FALSE
#> 8      8       -8348.43      TRUE         51  <NA>    FALSE
#> 9      9       -8348.43      TRUE         42  <NA>    FALSE
#> 10    10       -8348.43      TRUE         42  <NA>    FALSE
```

All ten starting values converged to the same solution, with no reported
boundary issues and a maximized log likelihood of approximately
-8348.429. This replication across starts provides evidence that the
fitted solution is numerically stable. The comparisons of alternative
profile and class counts below therefore begin from a solution that was
recovered consistently across multiple initializations.

## Classification uncertainty

An important aspect of model evaluation is the degree of uncertainty in
assigning observations to latent profiles or classes. Relative entropy
summarizes this uncertainty by measuring how concentrated the posterior
membership probabilities are. For $`N`$ units, $`K`$ latent classes, and
posterior probability $`p_{uk}`$ that unit $`u`$ belongs to class $`k`$,
relative entropy is defined as:

``` math
E=1-\frac{-\sum_{u=1}^{N}\sum_{k=1}^{K}p_{uk}\log p_{uk}}{N\log K}.
```

Relative entropy ranges from zero to one. A value of zero occurs when
each unit has equal posterior probability of belonging to every class,
which is maximal classification uncertainty. A value of one occurs when
every unit has posterior probability one for a single class and zero for
all others. Values between these extremes reflect varying degrees of
certainty in the posterior assignments.

In a multilevel latent profile model, relative entropy is calculated
separately for the enrolment and student levels. Enrolment-level entropy
summarizes the concentration of posterior probabilities for the
observation-level profiles, whereas student-level entropy summarizes
uncertainty in assignment to the higher-level latent classes. These
quantities should therefore be interpreted separately because
classification may be more distinct at one level than at the other.

Relative entropy is conditional on the fitted model and should not be
interpreted as a general measure of model correctness. A misspecified
model can still produce highly concentrated posterior probabilities. For
example, inappropriate assumptions about indicator independence or the
number of latent profiles may yield confident classifications even when
the model provides an inadequate representation of the data. Entropy
should therefore be considered alongside other diagnostic and
model-comparison measures. Relative entropy is also undefined at any
level containing only a single latent class because there is no
classification uncertainty to evaluate.

The command `get_results(fit, "entropy")` reports both the unnormalized
entropy sum and the corresponding relative entropy for the enrolment and
student levels. These results provide a descriptive summary of how
clearly the fitted model assigns units to the latent profiles and
classes at each level.

``` r

get_results(fit, "entropy")
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         3    1422    324.6527         0.792186
#> 2      groups         3     106     23.6028         0.797319
```

The fitted model has a relative entropy of 0.792 at the enrolment level
and 0.797 at the student level. These values indicate that posterior
membership probabilities are reasonably concentrated at both levels, but
that classification remains uncertain for some enrolments and students.
Consequently, reducing each posterior probability distribution to a
single most likely profile or class results in a loss of information.
Where possible, subsequent analyses should therefore retain the
uncertainty represented by the full posterior probabilities.

A complementary measure of classification certainty is the average
posterior probability within each assigned class. Let $`A_k`$ denote the
set of units whose most probable assignment is to class $`k`$. The
average posterior probability for that class is

``` math
\bar{p}_k
=
\frac{1}{|A_k|}
\sum_{u\in A_k} p_{uk},
```

where $`p_{uk}`$ is the posterior probability that unit $`u`$ belongs to
class $`k`$. This quantity summarizes, among units assigned to a
particular class, how strongly the model supports those assignments.
Values closer to one indicate that units assigned to the class generally
have high posterior probabilities for that class, whereas lower values
indicate greater ambiguity in the assignments.

Average posterior probabilities provide class-specific information that
is not visible from a single overall entropy measure. For example, a
model may show relatively high overall entropy while still containing
one profile or class whose members are considerably more difficult to
classify than those in the others. Examining average posterior
probabilities by class can therefore help identify where classification
uncertainty is concentrated.

Like relative entropy, however, average posterior probability is
conditional on the fitted model and its estimated classes. It describes
the certainty of the resulting modal assignments but does not establish
that the chosen number of classes, measurement assumptions, or
substantive interpretation is correct. Classification diagnostics should
therefore be considered alongside other evidence about model fit and
specification.

The command `get_results(fit, "classification")` reports these
class-specific classification summaries together with the number of
units assigned to each profile or class. The class counts are also
useful for identifying very small groups, for which classification
summaries and subsequent substantive comparisons may be less stable.

``` r

get_results(fit, "classification")
#>         level class n_modal proportion_modal estimated_n estimated_proportion average_posterior
#> 1 individuals     1     479         0.336850    478.0700             0.336196          0.852921
#> 2 individuals     2     576         0.405063    574.9629             0.404334          0.977847
#> 3 individuals     3     367         0.258087    368.9671             0.259471          0.845179
#> 4      groups     1      30         0.283019     31.0268             0.292706          0.977332
#> 5      groups     2      25         0.235849     27.6526             0.260874          0.897729
#> 6      groups     3      51         0.481132     47.3206             0.446421          0.864388
#>   odds_correct_classification
#> 1                    11.45003
#> 2                    65.02681
#> 3                    15.58023
#> 4                   104.18292
#> 5                    24.87018
#> 6                     7.90399
```

The class-specific results show that classification certainty varies
across profiles and student classes. At the enrolment level, 479
enrolments are assigned to the intermediate-activity Profile 1, with an
average posterior probability of 0.853. By comparison, the 576
enrolments assigned to the lower-activity Profile 2 have a much higher
average posterior probability of 0.978. This indicates that enrolments
classified into Profile 1 are, on average, more ambiguous than those
classified into Profile 2.

A similar pattern appears at the student level. Students assigned to
Class 3 have an average posterior probability of 0.864, whereas those
assigned to Class 1 have an average posterior probability of 0.977. The
student-level classifications are therefore not equally certain across
classes, even though the overall relative entropy is similar across the
two levels.

These differences matter for subsequent analyses because modal class
assignments treat the most likely class as if it were observed without
error. When average posterior probabilities are lower, that
simplification discards more of the uncertainty contained in the
posterior distribution. Later analyses should therefore retain or adjust
for classification uncertainty.

The [covariate
vignette](https://pak.dynasite.org/latents/articles/case-engagement-covariates.md)
describes approaches for incorporating this uncertainty when relating
latent profiles and classes to external variables, including methods
that adjust explicitly for classification error.

## Residual dependence

An important assumption of the diagonal covariance model is that the
indicators are conditionally independent within each latent profile.
Bivariate residual correlations can be used to examine whether
meaningful associations remain between pairs of indicators after the
fitted profile structure has been taken into account. For indicators
$`a`$ and $`b`$ within profile $`k`$, the residual correlation is
defined as

``` math
r_{ab,k}^{\text{observed}}
-
r_{ab,k}^{\text{expected}}.
```

The observed correlation is calculated using posterior probabilities of
profile membership as weights. Under the diagonal covariance model, the
expected within-profile correlation is zero. A residual close to zero
therefore indicates that the conditional-independence assumption
provides a reasonable approximation for that particular indicator pair.
Larger residuals indicate associations that remain unexplained by the
fitted measurement model, while the sign of the residual indicates
whether the remaining association is positive or negative.

These residual correlations are diagnostic quantities. Their values
depend on the estimated posterior profile memberships and therefore on
the fitted model itself. In addition, the reported p-values are
descriptive approximations. They do not account fully for uncertainty in
the estimated measurement parameters or posterior memberships, nor do
they adjust for the dependence created by repeated enrolments from the
same student. The magnitude and pattern of the residuals are therefore
more useful for identifying measurement assumptions that warrant further
examination than the p-values are for making formal inferential
decisions.

The command `get_results(fit, "residuals")` returns the indicator pairs
ordered by the absolute magnitude of their residual correlations.
Examining the largest residuals provides a convenient way to identify
the pairs for which the diagonal covariance assumption is least
consistent with the observed associations.

``` r

get_results(fit, "residuals") |> head(4)
#>     profile indicator_1 indicator_2     kind observed expected residual effective_n statistic df
#> 1 profile_2  forum_read  attendance gaussian 0.211091        0 0.211091     574.963   5.12547 NA
#> 2 profile_2      browse  attendance gaussian 0.177099        0 0.177099     574.963   4.28058 NA
#> 3 profile_2    lectures  attendance gaussian 0.113807        0 0.113807     574.963   2.73363 NA
#> 4 profile_2  forum_post  attendance gaussian 0.105264        0 0.105264     574.963   2.52683 NA
#>       p_value  p_adjusted
#> 1 2.96802e-07 2.96802e-07
#> 2 1.86407e-05 1.86407e-05
#> 3 6.26406e-03 6.26406e-03
#> 4 1.15097e-02 1.15097e-02
```

The largest residual correlation is 0.211 between forum reading and
attendance in the lower-activity Profile 2. The next largest is 0.177
between browsing and attendance, also in Profile 2. Notably, all four of
the largest residuals involve the attendance indicator. This pattern
suggests that the three-profile model does not completely account for
the association between recorded online activity and the number of days
on which students are active.

The concentration of the largest residuals around attendance provides a
specific direction for further model evaluation. The results identify a
particular part of the measurement structure where the
conditional-independence assumption may be too restrictive. Alternative
specifications could therefore examine whether allowing selected
within-profile associations involving attendance changes the estimated
profiles or improves the representation of the engagement data.

## Class counts and covariance assumptions

Class enumeration evaluates alternative models with different numbers of
enrolment-level profiles and student-level latent classes. One common
criterion for comparing these models is the Bayesian information
criterion (BIC),

``` math
\mathrm{BIC}
=
-2\ell + q\log(N),
```

where $`\ell`$ is the maximized log likelihood, $`q`$ is the number of
freely estimated parameters, and $`N`$ is the sample size used in the
penalty term. Lower BIC values indicate a more favorable balance between
model fit and complexity within the set of candidate models being
compared. In the present analysis, `bic_groups` uses the number of
students as $`N`$ for every candidate, so the BIC values are calculated
on a consistent basis across models.

The comparison also examines whether assumptions about the covariance
structure affect the preferred numbers of profiles and student classes.
With `model = "VVI"`, each profile has its own indicator variances, but
residual covariances between indicators are fixed to zero. This
specification corresponds to the conditional-independence assumption
examined in the preceding section. By contrast, `model = "VVV"`
estimates a separate unrestricted covariance matrix for each profile, so
the indicators can remain correlated within profiles.

Both covariance structures retain Gaussian enrolment-level profiles and
use common measurement parameters across student classes. The main
difference is therefore the amount of within-profile dependence that the
measurement model is permitted to capture. Because unrestricted
covariance matrices require substantially more parameters,
full-covariance models can be more difficult to estimate. Numerical
instability is particularly likely when some profiles contain relatively
few observations or when several likelihood maxima provide similar fits.

The
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
function fits every requested combination of profile count,
student-class count, and covariance structure to the same enrolments and
indicators. The resulting candidate models can then be summarized with
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
which reports BIC together with numerical diagnostics such as
convergence status, active parameter bounds, and replication of the best
likelihood across starting values.

These diagnostics should be considered before interpreting differences
in BIC. Comparisons are most informative when the candidate models have
converged, no relevant parameter bounds are active, and the best
likelihood is reproduced across multiple starts. BIC then provides
evidence about relative fit within the specified candidate set, but it
does not by itself establish that the Gaussian measurement assumptions
are correct or that a unique number of latent profiles or classes
exists.

``` r

candidates <- enumerate_classes(course_engagement, vars, id = "student",
                                n_profiles = 2:3, n_group_classes = 2:3,
                                model = c("VVI", "VVV"), n_starts = 6,
                                max_iter = 2000, tol = 1e-10, seed = 1)
#> Warning: Profile 2: 18.6 Gaussian membership mass for 20 free means and variances. This is a
#> conservative diagnostic, not an identification test; inspect the Gaussian covariance diagnostics.
get_results(candidates, "candidates")
#>               candidate n_profiles n_group_classes model log_likelihood n_parameters     aic
#> 1 profiles2_groups2_VVI          2               2   VVI       -8439.82           23 16925.6
#> 2 profiles3_groups2_VVI          3               2   VVI       -8364.52           35 16799.0
#> 3 profiles2_groups3_VVI          2               3   VVI       -8421.06           25 16892.1
#> 4 profiles3_groups3_VVI          3               3   VVI       -8348.43           38 16772.9
#> 5 profiles2_groups2_VVV          2               2   VVV       -8313.11           43 16712.2
#> 6 profiles3_groups2_VVV          3               2   VVV       -8289.07           65 16708.1
#> 7 profiles2_groups3_VVV          2               3   VVV       -8295.38           45 16680.8
#> 8 profiles3_groups3_VVV          3               3   VVV       -8275.08           68 16686.2
#>       kic bic_groups bic_individual sabic_groups sabic_individual caic_groups caic_individual
#> 1 16951.6    16986.9        17046.6      16914.2          16973.5     17009.9         17069.6
#> 2 16837.0    16892.3        16983.1      16781.7          16872.0     16927.3         17018.1
#> 3 16920.1    16958.7        17023.6      16879.7          16944.2     16983.7         17048.6
#> 4 16813.9    16874.1        16972.7      16754.0          16852.0     16912.1         17010.7
#> 5 16758.2    16826.8        16938.4      16690.9          16801.8     16869.8         16981.4
#> 6 16776.1    16881.3        17050.0      16675.9          16843.5     16946.3         17115.0
#> 7 16728.8    16800.6        16917.5      16658.4          16774.5     16845.6         16962.5
#> 8 16757.2    16867.3        17043.8      16652.4          16827.8     16935.3         17111.8
#>   awe_groups awe_individual icl_groups icl_individual clc_groups clc_individual profile_entropy
#> 1    17292.2        17411.7    17116.0        17175.7    17008.7        17008.7        0.938173
#> 2    17823.7        18005.5    17555.5        17646.4    17392.3        17392.3        0.789625
#> 3    17306.3        17436.1    17114.7        17179.6    16998.1        16998.1        0.940203
#> 4    17855.8        18053.1    17564.6        17663.3    17387.4        17387.4        0.792186
#> 5    17346.2        17569.5    17016.7        17128.4    16816.2        16816.2        0.906299
#> 6    17593.3        17930.8    17095.2        17263.9    16792.0        16792.0        0.933156
#> 7    17359.5        17593.2    17014.7        17131.5    16804.8        16804.8        0.909291
#> 8    17957.7        18310.8    17436.6        17613.2    17119.5        17119.5        0.828142
#>   group_entropy avepp_min avepp_max avepp_mean avepp_below group_avepp_min group_avepp_max
#> 1      0.939643  0.981757  0.984935   0.983346           0        0.968985        0.992014
#> 2      0.943149  0.848969  0.978772   0.892726           0        0.977176        0.992210
#> 3      0.816454  0.981777  0.985771   0.983774           0        0.900351        0.962213
#> 4      0.797319  0.845179  0.977847   0.891982           0        0.864388        0.977332
#> 5      0.942161  0.973255  0.977752   0.975503           0        0.970873        0.993019
#> 6      0.943398  0.796376  0.975771   0.915213           1        0.979379        0.990214
#> 7      0.816138  0.974333  0.977839   0.976086           0        0.903070        0.959803
#> 8      0.823957  0.649407  0.973627   0.840806           1        0.878064        0.955367
#>   group_avepp_mean group_avepp_below converged boundary min_standardized_eigenvalue
#> 1         0.980500                 0      TRUE    FALSE                   0.3714788
#> 2         0.984693                 0      TRUE    FALSE                   0.1999758
#> 3         0.922505                 0      TRUE    FALSE                   0.3748293
#> 4         0.913150                 0      TRUE    FALSE                   0.2002897
#> 5         0.981946                 0      TRUE    FALSE                   0.2710424
#> 6         0.984796                 0      TRUE    FALSE                   0.0110647
#> 7         0.922502                 0      TRUE    FALSE                   0.2720518
#> 8         0.920423                 0      TRUE    FALSE                   0.1201323
#>   min_gaussian_mass n_best_replicated
#> 1          587.9042                 6
#> 2          382.5560                 6
#> 3          591.0872                 6
#> 4          368.9671                 6
#> 5          583.2209                 6
#> 6           18.6125                 1
#> 7          589.4862                 6
#> 8          116.2557                 5
#>                                                                                                                                                                                warnings
#> 1                                                                                                                                                                                      
#> 2                                                                                                                                                                                      
#> 3                                                                                                                                                                                      
#> 4                                                                                                                                                                                      
#> 5                                                                                                                                                                                      
#> 6 Profile 2: 18.6 Gaussian membership mass for 20 free means and variances. This is a conservative diagnostic, not an identification test; inspect the Gaussian covariance diagnostics.
#> 7                                                                                                                                                                                      
#> 8                                                                                                                                                                                      
#>   error
#> 1  <NA>
#> 2  <NA>
#> 3  <NA>
#> 4  <NA>
#> 5  <NA>
#> 6  <NA>
#> 7  <NA>
#> 8  <NA>
```

Under the diagonal covariance specification, the lowest BIC is obtained
for the model with three profiles and three student classes, with a BIC
of 16874.1. When residual covariances are allowed to vary freely within
profiles, however, the preferred profile count changes. The two-profile,
three-class full-covariance model has a BIC of 16800.6, the lowest value
among the models in the candidate grid, and its best likelihood is
reproduced across all six starting values. The three-profile,
three-class full-covariance model has a higher BIC of 16867.3.

Most candidate models also show stable likelihood replication across
starts. The exception is the three-profile, two-class full-covariance
model, whose best solution is recovered by only one of the six
initializations. This weaker replication suggests greater numerical
uncertainty for that candidate and illustrates why class enumeration
should not be based on information criteria alone.

Taken together, these results show that support for the original
three-profile, three-class specification depends on the assumptions made
about within-profile covariance. Under the diagonal model, the
three-by-three solution is preferred within the candidate set. Once
residual correlations are estimated directly, a two-profile, three-class
model receives stronger BIC support. The original three-by-three model
should therefore be regarded as a working representation whose support
depends on the measurement assumptions.

The
[`get_candidate()`](https://pak.dynasite.org/latents/reference/get_candidate.md)
function can be used to retrieve a fitted candidate with a particular
combination of profile count, student-class count, and covariance
structure. Examining its residual correlations provides a direct check
of whether the estimated covariance parameters reproduce the
within-profile associations observed in the data.

``` r

full <- get_candidate(candidates, n_profiles = 2, n_group_classes = 3,
                      model = "VVV")
get_results(full, "residuals") |> head(4)
#>     profile indicator_1 indicator_2     kind   observed  expected    residual effective_n
#> 1 profile_1      browse  forum_read gaussian 0.04220405 0.0422031 8.99773e-07     832.514
#> 2 profile_1  forum_read  attendance gaussian 0.22563970 0.2256388 8.75914e-07     832.514
#> 3 profile_1    lectures  forum_read gaussian 0.00953707 0.0095362 8.64190e-07     832.514
#> 4 profile_2  forum_read  forum_post gaussian 0.07712515 0.0771244 7.70121e-07     589.486
#>     statistic df  p_value p_adjusted
#> 1 2.59609e-05 NA 0.999979   0.999979
#> 2 2.65807e-05 NA 0.999979   0.999979
#> 3 2.48920e-05 NA 0.999980   0.999980
#> 4 1.87620e-05 NA 0.999985   0.999985
```

In the two-profile, three-class full-covariance model, the estimated
correlation between forum reading and attendance in Profile 1 is 0.226,
which matches the corresponding observed within-profile correlation to
six decimal places. As a result, even the largest residual correlations
are close to zero. This is expected because the full-covariance
specification estimates these within-profile associations directly
rather than constraining them to zero.

The close agreement between the observed and expected correlations
therefore verifies that the fitted covariance structure is reproducing
the associations it was designed to model. It should not, however, be
interpreted as independent evidence that the two-profile solution or the
three student classes are substantively correct. Those questions require
consideration of the broader pattern of model fit, numerical stability,
classification quality, and substantive interpretability.

## Uncertainty in profile means

Parameter inference quantifies the sampling uncertainty associated with
an estimated model parameter, $`\hat{\theta}`$. By default, the standard
error is obtained from the observed information matrix. Specifically,
the estimated variance of a parameter is taken from the corresponding
diagonal element of the inverse observed information matrix, and the
standard error is the square root of that value.

A conventional 95% Wald confidence interval is calculated as

``` math
\hat{\theta}
\pm
1.96\,\operatorname{SE}(\hat{\theta}).
```

These intervals provide a local approximation to parameter uncertainty
within the fitted model. Their use requires a satisfactorily converged
solution, inactive parameter bounds, and an invertible information
matrix. If these conditions are not met, the resulting standard errors
and confidence intervals may be unstable or unavailable.

The intervals are also conditional on the specified latent structure. In
particular, they treat the selected numbers of enrolment-level profiles
and student-level classes, together with the assumed covariance
structure, as fixed. A parameter estimate can therefore have a narrow
confidence interval even when an alternative model specification
supports a different number of profiles. Wald intervals describe
uncertainty about parameters within a chosen model; they do not
incorporate uncertainty about the choice of model itself.

The
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
function calculates standard errors and confidence intervals for every
parameter of the original three-profile, three-class model. The
discussion below uses the rows for the attendance means.

``` r

parameter_inference(fit)
#>          level       outcome          term   parameter   estimate standard_error statistic
#> 1  measurement     profile_1        browse        mean  0.3870324      0.0454157   8.52199
#> 2  measurement     profile_1      lectures        mean  0.1856519      0.0575759   3.22447
#> 3  measurement     profile_1    forum_read        mean  0.4464107      0.0441603  10.10886
#> 4  measurement     profile_1    forum_post        mean  0.3376538      0.0533639   6.32738
#> 5  measurement     profile_1    attendance        mean  0.2564040      0.0515939   4.96966
#> 6  measurement     profile_2        browse        mean -0.7768095      0.0314788 -24.67724
#> 7  measurement     profile_2      lectures        mean -0.6123273      0.0312717 -19.58085
#> 8  measurement     profile_2    forum_read        mean -0.8914204      0.0264490 -33.70342
#> 9  measurement     profile_2    forum_post        mean -0.7623134      0.0280090 -27.21677
#> 10 measurement     profile_2    attendance        mean -0.9399557      0.0268401 -35.02059
#> 11 measurement     profile_3        browse        mean  0.7089196      0.0514602  13.77606
#> 12 measurement     profile_3      lectures        mean  0.7134804      0.0616461  11.57381
#> 13 measurement     profile_3    forum_read        mean  0.8109612      0.0462641  17.52897
#> 14 measurement     profile_3    forum_post        mean  0.7503371      0.0525223  14.28607
#> 15 measurement     profile_3    attendance        mean  1.1328124      0.0610244  18.56326
#> 16 measurement     profile_1        browse    variance  0.6132450      0.0450473        NA
#> 17 measurement     profile_1      lectures    variance  0.7610659      0.0586588        NA
#> 18 measurement     profile_1    forum_read    variance  0.4571729      0.0359665        NA
#> 19 measurement     profile_1    forum_post    variance  0.7116451      0.0561877        NA
#> 20 measurement     profile_1    attendance    variance  0.1957470      0.0216178        NA
#> 21 measurement     profile_2        browse    variance  0.5240143      0.0322067        NA
#> 22 measurement     profile_2      lectures    variance  0.5359582      0.0321819        NA
#> 23 measurement     profile_2    forum_read    variance  0.3571087      0.0226928        NA
#> 24 measurement     profile_2    forum_post    variance  0.4139024      0.0262015        NA
#> 25 measurement     profile_2    attendance    variance  0.3626838      0.0227638        NA
#> 26 measurement     profile_3        browse    variance  0.5195008      0.0482432        NA
#> 27 measurement     profile_3      lectures    variance  0.8071150      0.0690787        NA
#> 28 measurement     profile_3    forum_read    variance  0.4642238      0.0394977        NA
#> 29 measurement     profile_3    forum_post    variance  0.5837024      0.0553610        NA
#> 30 measurement     profile_3    attendance    variance  0.2025563      0.0246588        NA
#> 31     profile     profile_1 group_class_1 probability  0.0701381      0.0205420        NA
#> 32     profile     profile_2 group_class_1 probability  0.8466733      0.0245137        NA
#> 33     profile     profile_3 group_class_1 probability  0.0831886      0.0188275        NA
#> 34     profile     profile_1 group_class_2 probability  0.5470108      0.0746365        NA
#> 35     profile     profile_2 group_class_2 probability  0.0552359      0.0247742        NA
#> 36     profile     profile_3 group_class_2 probability  0.3977533      0.0757776        NA
#> 37     profile     profile_1 group_class_3 probability  0.3881406      0.0428633        NA
#> 38     profile     profile_2 group_class_3 probability  0.3171672      0.0364866        NA
#> 39     profile     profile_3 group_class_3 probability  0.2946922      0.0407032        NA
#> 40       group group_class_1          <NA> probability  0.2927056      0.0482765        NA
#> 41       group group_class_2          <NA> probability  0.2608736      0.0741567        NA
#> 42       group group_class_3          <NA> probability  0.4464207      0.0737344        NA
#>         p_value   p_adjusted   conf_low conf_high
#> 1   1.56830e-17  1.56830e-17  0.2980192  0.476046
#> 2   1.26205e-03  1.26205e-03  0.0728052  0.298499
#> 3   5.04657e-24  5.04657e-24  0.3598581  0.532963
#> 4   2.49357e-10  2.49357e-10  0.2330625  0.442245
#> 5   6.70708e-07  6.70708e-07  0.1552819  0.357526
#> 6  1.87769e-134 1.87769e-134 -0.8385067 -0.715112
#> 7   2.25242e-85  2.25242e-85 -0.6736188 -0.551036
#> 8  5.15138e-249 5.15138e-249 -0.9432594 -0.839581
#> 9  4.11293e-163 4.11293e-163 -0.8172099 -0.707417
#> 10 1.09354e-268 1.09354e-268 -0.9925613 -0.887350
#> 11  3.55114e-43  3.55114e-43  0.6080593  0.809780
#> 12  5.59425e-31  5.59425e-31  0.5926563  0.834305
#> 13  8.61179e-69  8.61179e-69  0.7202853  0.901637
#> 14  2.67283e-46  2.67283e-46  0.6473953  0.853279
#> 15  6.37199e-77  6.37199e-77  1.0132067  1.252418
#> 16           NA           NA  0.5310154  0.708208
#> 17           NA           NA  0.6543593  0.885173
#> 18           NA           NA  0.3918458  0.533391
#> 19           NA           NA  0.6096171  0.830749
#> 20           NA           NA  0.1576488  0.243052
#> 21           NA           NA  0.4645441  0.591098
#> 22           NA           NA  0.4764530  0.602895
#> 23           NA           NA  0.3152899  0.404474
#> 24           NA           NA  0.3656065  0.468578
#> 25           NA           NA  0.3207027  0.410160
#> 26           NA           NA  0.4330518  0.623208
#> 27           NA           NA  0.6824699  0.954525
#> 28           NA           NA  0.3929202  0.548467
#> 29           NA           NA  0.4846849  0.702948
#> 30           NA           NA  0.1595593  0.257140
#> 31           NA           NA  0.0390943  0.122686
#> 32           NA           NA  0.7922662  0.888830
#> 33           NA           NA  0.0529689  0.128313
#> 34           NA           NA  0.4008901  0.685456
#> 35           NA           NA  0.0225372  0.129110
#> 36           NA           NA  0.2621452  0.551115
#> 37           NA           NA  0.3081288  0.474675
#> 38           NA           NA  0.2502982  0.392547
#> 39           NA           NA  0.2215729  0.380157
#> 40           NA           NA  0.2076223  0.395263
#> 41           NA           NA  0.1424337  0.428581
#> 42           NA           NA  0.3100415  0.591372
```

For the higher-activity Profile 3, the estimated attendance mean is
1.133 standard deviations above the course mean, with a standard error
of 0.061 and a 95% Wald interval from 1.013 to 1.252. The
intermediate-activity Profile 1 has an interval from 0.155 to 0.358,
while the lower-activity Profile 2 has an interval from $`-0.993`$ to
$`-0.887`$.

Within the original three-profile, three-class specification, these
intervals show clear separation among the attendance means. This
provides evidence that the fitted profiles differ substantially in their
expected levels of attendance under that model. However, the precision
of these estimates should not be interpreted as resolving the
model-specification uncertainty identified in the covariance comparison.
The confidence intervals quantify uncertainty in the profile means
conditional on the three-by-three model; they do not address whether a
different covariance structure or profile count provides a preferable
representation of the data.

## Limits of the interpretation

The model comparisons presented here should be interpreted within the
limits of the candidate set that was examined. Only a restricted range
of profile counts, student-class counts, and covariance structures was
considered, and the resulting BIC values also depend on the convention
used to define the sample size in the BIC penalty. The comparison
therefore identifies relative support among the fitted candidates only.

Additional caution is warranted for candidate models whose best
likelihood was not reproduced across multiple starting values. In
particular, full-covariance solutions with weak replication should be
refitted with additional starts before they are used for substantive
interpretation. More generally, all reported classification summaries,
residual diagnostics, standard errors, and confidence intervals are
conditional on the measurement specification from which they were
obtained.

The multilevel structure accounts for the nesting of repeated enrolments
within students by allowing students to differ in the composition of
their engagement profiles. However, the model does not directly estimate
serial dependence between adjacent enrolments. The argument
`time = "sequence"` preserves the ordering of each student’s enrolments
and makes that ordering available for sequence-based summaries, but it
does not by itself introduce a transition or autoregressive structure
into the fitted model.

When the scientific question concerns how students move from one
engagement profile to another over time, a transition model is more
appropriate. The [transition
vignette](https://pak.dynasite.org/latents/articles/case-engagement-transitions.md)
extends the analysis by estimating probabilities of movement between
profiles across successive enrolments.
