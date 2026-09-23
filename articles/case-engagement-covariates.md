# Case study: Covariates and staged estimation in multilevel LPA

Three-step covariate analysis is a method for examining associations
between latent profile membership and variables that are external to the
measurement model. An important feature of this approach is that it
adjusts for uncertainty in profile classification. The adjustment
reduces the bias that arises when individuals are assigned to latent
classes with imperfect certainty.

In the course-engagement example, the primary research question is
whether students’ previous grades are associated with their engagement
patterns in a subsequent course. This relationship can be investigated
in several ways. One approach is to regress latent profile membership on
previous grades in order to assess whether prior academic performance
predicts the likelihood of belonging to a particular engagement profile.
A second approach is to compare mean previous grades across the
identified profiles. A third approach is to estimate the model jointly
by including previous grades directly in the latent profile model.

These approaches differ in the extent to which the external variable is
allowed to influence the latent profile solution. In the three-step
approaches, the measurement model is estimated first, and the resulting
profile definitions are held fixed while associations with previous
grades are examined. This preserves the original interpretation of the
engagement profiles. In one-step estimation, previous grade enters the
latent profile model and can change the profiles. The final section
turns to staged estimation, which holds the profile definitions fixed
while the student classes are estimated.

``` r

library(latents)
```

## Data and engagement profiles

The bundled synthetic `course_engagement` dataset contains enrolment
records nested within students. Engagement is represented by five
activity indicators: page views (`browse`), lecture videos viewed
(`lectures`), forum posts read (`forum_read`), forum posts written
(`forum_post`), and active days (`attendance`). Each indicator is
simulated on the `log1p` scale of an activity count and standardized
within course to make activity levels comparable across courses.

The variable `previous_grade` represents the standardized grade obtained
in the preceding course. For a student’s first observed course, an entry
grade is used instead. Because this grade is measured before engagement
in the subsequent course, it can be treated as a predictor of later
engagement patterns. The [introductory
vignette](https://pak.dynasite.org/latents/articles/case-engagement-profiles.md)
provides additional explanation of the longitudinal structure of the
data and the distinction between the student and enrolment levels.

The indicator vector `vars` specifies the activity measures used to
estimate the latent engagement profiles. In keeping `previous_grade`
separate from these indicators, the engagement profiles can first be
established using activity data alone. The relationship between previous
academic performance and profile membership can then be examined in a
subsequent stage without allowing grades to influence the initial
definition of the profiles.

``` r

vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")
```

The
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
function estimates the means and variances that define each latent
profile, together with the probability of observing each profile within
each student class. Let this conditional probability be written as

``` math
\pi_{k\mid m}=P(C_{it}=k\mid G_i=m),
```

where $`C_{it}`$ denotes the latent engagement profile of student $`i`$
in course $`t`$, and $`G_i`$ denotes the student’s latent class. For
each student class, the profile probabilities sum to one across all
profiles. These probabilities therefore describe how the distribution of
course-level engagement profiles differs across student-level classes.

The model considered here includes three engagement profiles and three
student classes. Within each profile, the five engagement indicators are
modeled as independent Gaussian variables. The profile definitions are
held common across student classes, meaning that the same
profile-specific means and variances apply to all students. Student
classes differ instead in the probabilities with which the engagement
profiles occur.

Model estimation relies on several assumptions. Students are assumed to
be independent of one another, while enrolments are assumed to be
conditionally independent given student class. The model also assumes
conditional independence among the engagement indicators within each
profile. Violations of this latter assumption can affect the profiles
that are recovered from the data. In addition, because latent profile
models use likelihood-based estimation, the optimization procedure may
converge to different local maxima depending on its starting values.

The [evaluation
vignette](https://pak.dynasite.org/latents/articles/case-engagement-evaluation.md)
discusses indicator dependence, likelihood maxima, and related
model-evaluation issues in greater detail. In the present tutorial, the
same candidate model is retained throughout so that the different
covariate methods can be compared on a common basis. Ten random starts
are used to reduce sensitivity to local solutions, and a fixed random
seed ensures that the numerical results are reproducible.

``` r

fit <- multilpa(course_engagement, vars, id = "student", time = "sequence",
                n_profiles = 3, n_group_classes = 3, n_starts = 10,
                max_iter = 2000, tol = 1e-10, seed = 1)
fit
#> Two-level latent profile analysis: 3 profiles, 3 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -8348.428714 | AIC: 16772.857 | BIC (groups): 16874.068
#> Converged: TRUE | iterations: 164 | best start: 4/10
#> 
#>  profile  browse lectures forum_read forum_post attendance count proportion
#>        1  0.7088   0.7134     0.8109     0.7502     1.1326 369.2     0.2596
#>        2  0.3870   0.1855     0.4463     0.3375     0.2562 477.9     0.3361
#>        3 -0.7768  -0.6123    -0.8914    -0.7623    -0.9400 575.0     0.4043
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

The three profiles can be characterized as representing higher,
intermediate, and lower levels of course activity, respectively. This
ordering is particularly evident in the estimated means for
`attendance`, which are 1.133 for Profile 1, 0.256 for Profile 2, and
-0.940 for Profile 3. On this basis, Profile 3 is treated as the
lower-activity reference category when interpreting the
profile-membership regressions presented below.

These labels describe the current fitted solution; the profile numbers
themselves carry no fixed meaning. Because latent profile labels are
determined by the estimated parameters in each model, the
profile-specific means should be examined again whenever a new model is
fitted to confirm that the same substantive interpretation still
applies.

## Accounting for classification error

The three-step approach accounts for uncertainty in profile assignment
by using a classification-error matrix. This matrix describes the
probability that an individual is assigned to a particular profile given
their underlying latent profile. Let

``` math
D_{ks}=P(\widehat{C}=s\mid C=k),
```

where $`C`$ denotes the underlying latent profile and $`\widehat{C}`$
denotes the assigned profile, defined as the profile with the largest
posterior probability. Each row of $`D`$ sums to one because it
represents the distribution of possible assignments for a given
underlying profile. The diagonal elements therefore give the estimated
probabilities of being assigned to the correct profile, whereas the
off-diagonal elements represent misclassification between profiles. In
`multilpa`, this matrix can be extracted with
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
using `"classification_errors"`.

The correction relies on the assumption that, conditional on the
underlying latent profile, classification errors are independent of the
auxiliary variable being studied. In other words, once true profile
membership is taken into account, the probability of being misclassified
is assumed not to depend systematically on the covariate or distal
outcome.

The classification-error matrix is estimated from the initial
measurement model and is subsequently treated as fixed during the
covariate analysis. As a result, confidence intervals and standard
errors from the later stages do not incorporate uncertainty associated
with estimating this matrix. This distinction is important when
interpreting inferential results, particularly when profile separation
is limited. Poorly separated profiles produce greater classification
uncertainty and can make the error correction less stable.

``` r

get_results(fit, "classification_errors", level = "individuals")
#>         level true_class assigned_class probability
#> 1 individuals          1              1   8.418e-01
#> 2 individuals          1              2   1.580e-01
#> 3 individuals          1              3   1.515e-04
#> 4 individuals          2              1   1.197e-01
#> 5 individuals          2              2   8.537e-01
#> 6 individuals          2              3   2.660e-02
#> 7 individuals          3              1   5.238e-06
#> 8 individuals          3              2   2.038e-02
#> 9 individuals          3              3   9.796e-01
```

The estimated probabilities of correct profile assignment are 0.842 for
Profile 1, 0.854 for Profile 2, and 0.980 for Profile 3. Classification
is therefore most precise for Profile 3, whereas greater uncertainty
remains between Profiles 1 and 2.

The main source of ambiguity is the overlap between these first two
profiles. The estimated probability of assigning an observation from
Profile 1 to Profile 2 is 0.158, while the probability of assigning an
observation from Profile 2 to Profile 1 is 0.120. These off-diagonal
probabilities indicate that a non-negligible proportion of observations
may be classified into the neighboring profile rather than their
underlying profile.

Treating the assigned profiles as if they were observed without error
would ignore this uncertainty. The three-step correction instead carries
the estimated misclassification structure forward into the subsequent
covariate analysis, so associations with external variables are
estimated with classification error taken into account.

## Previous grade as a membership predictor

The [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
function estimates the association between an external predictor and
latent profile membership while keeping the profile definitions and
classification-error matrix fixed. In this example, `previous_grade` is
used to predict membership in the engagement profiles identified in the
measurement model. Because Profile 3 represents the lower-activity
profile, it is used as the reference category.

Let $`x`$ denote standardized previous grade. For Profiles 1 and 2, the
multinomial regression can be written as

``` math
\log\frac{P(C=k\mid x)}{P(C=3\mid x)}
=
\alpha_k+\beta_k x,
\qquad k=1,2.
```

The intercept $`\alpha_k`$ represents the log odds of belonging to
Profile $`k`$, rather than the lower-activity reference profile, when
$`x=0`$. Because `previous_grade` is standardized, $`x=0`$ corresponds
to the average previous grade. The slope $`\beta_k`$ represents the
change in these log odds associated with a one-standard-deviation
increase in previous grade. A positive slope therefore indicates that
higher previous grades are associated with greater relative odds of
belonging to Profile $`k`$ rather than the lower-activity profile.

The model assumes that previous grade has a linear effect on the
log-odds scale. For analyses conducted at the enrolment level,
`vcov_type = "robust"` can be used to account for the clustering of
multiple enrolments within the same student. Under this option, score
contributions are clustered by student when estimating the covariance
matrix. Reliable estimation therefore requires a sufficient number of
independent students.

The [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
regression describes the marginal association between previous grade and
profile membership. This differs from the one-step model considered
later, in which profile membership is modeled conditional on student
class. The two approaches therefore answer related but distinct
questions about how previous academic performance is associated with
subsequent engagement.

The reported 95% confidence intervals are based on the normal
approximation. When several slope coefficients are tested
simultaneously, the reported adjusted $`p`$-values use the
Benjamini–Hochberg procedure to control the false discovery rate across
those tests. The confidence intervals remain pointwise intervals and are
not adjusted for multiple comparisons.

``` r

r3step(fit, data = course_engagement, covariates = "previous_grade", vcov_type = "robust")
#>         level outcome           term estimate standard_error statistic   p_value p_value_adjusted
#> 1 individuals class_1    (Intercept)  -0.4456        0.13637   -3.2676 1.085e-03               NA
#> 2 individuals class_1 previous_grade   0.6769        0.08623    7.8506 4.142e-15        8.284e-15
#> 3 individuals class_2    (Intercept)  -0.1416        0.14708   -0.9626 3.358e-01               NA
#> 4 individuals class_2 previous_grade   0.4753        0.06903    6.8854 5.762e-12        5.762e-12
#>   conf_low conf_high
#> 1  -0.7129   -0.1783
#> 2   0.5079    0.8459
#> 3  -0.4298    0.1467
#> 4   0.3400    0.6106
```

In the fitted model, a one-standard-deviation increase in previous grade
is associated with a 0.677 increase in the log odds of belonging to the
higher-activity profile rather than the lower-activity profile. The
corresponding 95% confidence interval ranges from 0.508 to 0.846. For
the intermediate-activity profile, the estimated increase in log odds is
0.475, with a 95% confidence interval from 0.340 to 0.611.

Both coefficients are positive: higher previous grades are associated
with greater relative odds of membership in the higher- or
intermediate-activity profiles compared with the lower-activity profile.
The results therefore show a consistent association between stronger
previous academic performance and more active engagement in the
subsequent course.

These findings are observational associations. Although `previous_grade`
is measured before engagement in the subsequent course, temporal
ordering alone does not establish that higher grades cause students to
adopt a more active engagement pattern.

## Comparing previous grades across profiles

The
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
function estimates profile-specific means of a numeric auxiliary
variable and tests differences between those means while accounting for
uncertainty in profile classification. In this example, `previous_grade`
is treated as the auxiliary variable, so prior academic performance can
be compared across the engagement profiles identified in the measurement
model.

By default,
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
uses the Bolck–Croon–Hagenaars (BCH) method. This approach constructs
signed weights from the inverse of the classification-error matrix and
uses those weights to correct for imperfect profile assignment. Let
$`x_u`$ denote the previous grade for enrolment $`u`$, and let
$`w_{uk}`$ denote its BCH weight for Profile $`k`$. The estimated
profile-specific mean is

``` math
\widehat{\mu}_k
=
\frac{\sum_u w_{uk}x_u}
{\sum_u w_{uk}}.
```

Differences between profiles are then expressed as pairwise contrasts.
For Profiles $`k`$ and $`l`$, the contrast is

``` math
\widehat{\mu}_k-\widehat{\mu}_l.
```

The BCH correction requires the classification-error matrix to be
invertible and the resulting weighted profile sizes to be sufficiently
well behaved for estimation. The default standard errors account for the
clustering of multiple enrolments within students, while the previously
estimated classification-error matrix is treated as fixed. Consequently,
uncertainty from estimation of that matrix is not propagated into the
reported standard errors.

Although
[`three_step()`](https://pak.dynasite.org/latents/reference/three_step.md)
can also be used with distal outcomes measured after profile membership,
`previous_grade` occurs before the engagement measures in this example.
The estimated contrasts should therefore be interpreted as differences
in prior academic performance across subsequently observed engagement
profiles rather than as effects of profile membership on grades.

The reported 95% confidence intervals are pointwise normal intervals.
Adjusted $`p`$-values use the Benjamini–Hochberg procedure across the
set of pairwise profile contrasts.

``` r

three_step(
  fit,
  data = course_engagement,
  outcome = "previous_grade",
  contrast = "pairs"
)
#>         level method class reference_class estimate standard_error statistic   p_value
#> 1 individuals    bch     2               1  -0.1828        0.09448    -1.935 5.304e-02
#> 2 individuals    bch     3               1  -0.6287        0.07750    -8.112 4.960e-16
#> 3 individuals    bch     3               2  -0.4459        0.06048    -7.374 1.661e-13
#>   p_value_adjusted conf_low conf_high
#> 1        5.304e-02  -0.3680  0.002397
#> 2        1.488e-15  -0.7806 -0.476828
#> 3        2.491e-13  -0.5645 -0.327412
```

The BCH-adjusted estimates show that previous grades in the
lower-activity profile are 0.629 standard deviations below those in the
higher-activity profile and 0.446 standard deviations below those in the
intermediate-activity profile. These contrasts indicate that students
who later exhibit lower engagement tend to have had lower grades in the
preceding course.

The difference between the intermediate- and higher-activity profiles is
smaller. The estimated Profile 2 versus Profile 1 contrast is
$`-0.183`$, with a 95% confidence interval from $`-0.368`$ to $`0.002`$.
Because this interval extends slightly across zero, the difference in
previous grades between the two more active profiles is less precisely
estimated than their respective differences from the lower-activity
profile.

### Comparison with modal assignment

A useful comparison is provided by modal analysis, in which each
enrolment is assigned with certainty to the profile having the largest
posterior probability. The same profile means and pairwise contrasts are
then calculated using these hard assignments rather than
classification-error-adjusted weights.

``` r

three_step(
  fit,
  data = course_engagement,
  outcome = "previous_grade",
  contrast = "pairs",
  method = "modal"
)
#>         level method class reference_class estimate standard_error statistic   p_value
#> 1 individuals  modal     2               1  -0.1430        0.06736    -2.123 3.379e-02
#> 2 individuals  modal     3               1  -0.5904        0.06669    -8.853 8.551e-19
#> 3 individuals  modal     3               2  -0.4474        0.05201    -8.602 7.833e-18
#>   p_value_adjusted conf_low conf_high
#> 1        3.379e-02  -0.2750  -0.01095
#> 2        2.565e-18  -0.7211  -0.45969
#> 3        1.175e-17  -0.5494  -0.34548
```

Modal analysis is simple to interpret, but it treats estimated profile
assignments as if they were observed without error. Its standard errors
therefore do not incorporate uncertainty associated with imperfect
classification.

The effect of this simplification is particularly visible for Profiles 1
and 2, which showed the greatest classification overlap. Under modal
assignment, the estimated Profile 2 versus Profile 1 difference is
$`-0.143`$, compared with $`-0.183`$ under the BCH correction. The
corresponding standard error is 0.067 under modal assignment and 0.094
under BCH.

Thus, accounting for classification error affects both the estimated
contrast and its estimated uncertainty. The difference is most
noticeable for the two profiles that are hardest to distinguish.
Classification adjustment matters most when external variables are
compared across partially overlapping latent profiles.

## Estimating profiles and grade associations together

One-step estimation differs from the staged approaches considered above
because the latent profile measurement model and the regression of
profile membership on the covariate are estimated simultaneously. In
this example, specifying `profile_covariates = "previous_grade"` in
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
allows information about previous grades to contribute directly to
estimation of the engagement profiles as well as to the association
between grades and profile membership.

Let $`x`$ denote standardized previous grade and $`G=m`$ denote
membership in student class $`m`$. With Profile 3 as the reference
category, the profile-membership model is

``` math
\log\frac{P(C=k\mid G=m,x)}
{P(C=3\mid G=m,x)}
=
\alpha_{km}+\beta_k x,
\qquad k=1,2.
```

The intercepts $`\alpha_{km}`$ are allowed to vary across student
classes. They therefore describe differences in the baseline profile
probabilities between student classes when previous grade is at its
standardized mean, $`x=0`$. In contrast, the grade slopes $`\beta_k`$
are shared across student classes. Each slope represents the change in
the log odds of belonging to Profile $`k`$, rather than Profile 3,
associated with a one-standard-deviation increase in previous grade,
conditional on student class.

The one-step model retains the same Gaussian measurement assumptions
used in the original latent profile model, including conditional
independence of the engagement indicators within profiles. The important
difference is that previous grade is now part of the model during
estimation rather than being introduced only after the profiles have
been established.

As a consequence, the profiles obtained from the one-step model need not
be identical to those from the original measurement model. Allowing
previous grade to contribute to estimation can alter profile-specific
means, variances, or membership probabilities. The estimated regression
coefficients also have a different interpretation from the three-step
coefficients because they describe the association between previous
grade and profile membership conditional on student class, whereas the
earlier
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md)
analysis describes marginal profile membership.

For these reasons, the magnitude of the one-step coefficients should not
be expected to match the corresponding three-step estimates exactly.
Before interpreting the coefficients, the fitted profile means should be
inspected to determine whether Profiles 1, 2, and 3 still correspond to
higher-, intermediate-, and lower-activity engagement patterns. Profile
numbers are labels generated by the fitted model and do not carry fixed
substantive meanings across separate model fits.

``` r

covariate_fit <- multilpa(
  course_engagement,
  vars,
  id = "student",
  n_profiles = 3,
  n_group_classes = 3,
  profile_covariates = "previous_grade",
  n_starts = 6,
  max_iter = 2000,
  tol = 1e-10,
  seed = 1
)
#> Warning: The best covariate start had not converged when `max_iter` was reached, or its membership
#> logits still carry a non-negligible score, so the returned estimate is not a maximum.

covariate_fit
#> Multilevel LPA with covariates: 3 profiles, 3 group classes
#> Log likelihood -8331.180650; AIC 16742.361; BIC (groups) 16848.899; converged FALSE
#> 
#>  profile  browse lectures forum_read forum_post attendance count proportion
#>        1  0.7278   0.7277     0.8141     0.7569     1.1587 347.3     0.2442
#>        2  0.3881   0.2029     0.4610     0.3536     0.2772 498.4     0.3505
#>        3 -0.7743  -0.6141    -0.8891    -0.7620    -0.9379 576.3     0.4053
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

This fit reports `converged FALSE`, and
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
warns that the best start had not converged when `max_iter` was reached,
so the returned estimate is not a maximum. The profile means and
membership coefficients of a one-step model are meaningful only when the
fit has converged. They are therefore not interpreted here; a
substantive analysis would first resolve the convergence problem, for
example with more starts or iterations.

The `"coefficients"` table shows the form of the result: one intercept
for each profile and student class, one `previous_grade` slope for each
profile relative to the reference profile, and the intercepts of the
student-class model.

``` r

get_results(covariate_fit, "coefficients")
#>      level       outcome           term parameter estimate
#> 1  profile     profile_1  group_class_1     logit -2.25650
#> 2  profile     profile_1  group_class_2     logit -0.16147
#> 3  profile     profile_1  group_class_3     logit  1.68877
#> 4  profile     profile_1 previous_grade     logit  0.53323
#> 5  profile     profile_2  group_class_1     logit -2.52746
#> 6  profile     profile_2  group_class_2     logit  0.22908
#> 7  profile     profile_2  group_class_3     logit  2.25912
#> 8  profile     profile_2 previous_grade     logit  0.29421
#> 9    group group_class_1    (Intercept)     logit  0.07867
#> 10   group group_class_2    (Intercept)     logit  0.53863
```

In a converged fit, a positive slope $`\beta_k`$ indicates higher odds
of profile $`k`$ relative to profile 3 for each standard deviation of
previous grade, conditional on student class. Because these slopes are
conditional on student class and the profiles are estimated together
with the covariate, they need not equal the three-step slopes even when
both fits have converged.

## Holding profile definitions fixed across student classes

Staged estimation separates the estimation of the course-level
engagement profiles from the estimation of the higher-level
student-class structure. In the first stage, a common set of engagement
profiles is estimated without allowing student classes to alter the
profile measurement parameters. In the second stage, the composition of
student classes is estimated while the profile means and variances from
the first stage are held fixed.

The
[`fit_staged()`](https://pak.dynasite.org/latents/reference/fit_staged.md)
function implements this procedure by first fitting a model with a
single student class and then estimating the requested three-class
student-level structure. During the second stage, the model estimates
the profile probabilities within each student class,

``` math
\pi_{k\mid m}=P(C_{it}=k\mid G_i=m),
```

together with the marginal student-class probabilities,

``` math
\omega_m=P(G_i=m).
```

The profile-specific means and variances remain fixed at their
first-stage estimates. Staged estimation therefore isolates the
estimation of the higher-level class structure from the estimation of
the profile measurement model.

This separation has an important inferential consequence. Uncertainty in
the first-stage measurement parameters is not propagated into the
second-stage estimates. The resulting model is therefore conditional on
the fixed first-stage profile solution and should not be interpreted as
the joint maximum-likelihood solution that would be obtained by
estimating all parameters simultaneously.

Profile numbering also requires careful attention when comparing staged
and joint models. Latent profile and student-class labels are arbitrary
numerical identifiers, and separate estimation procedures may assign
different numbers to substantively similar groups. Comparisons should
therefore be based on the estimated activity patterns and class
compositions rather than on the numerical labels alone.

``` r

staged <- fit_staged(
  course_engagement,
  vars,
  id = "student",
  n_profiles = 3,
  n_group_classes = 3,
  n_starts = 6,
  max_iter = 2000,
  tol = 1e-10,
  seed = 1
)

staged
#> Two-level latent profile analysis: 3 profiles, 3 group classes
#> 1422 individuals in 106 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -8351.807729 | AIC: 16719.615 | BIC (groups): 16740.923
#> Converged: TRUE | iterations: 58 | best start: 2/6
#> Held fixed, not estimated here: means, variances (first stage)
#> Parameters estimated here: 8; with the held measurement: 38
#> 
#>  profile  browse lectures forum_read forum_post attendance count proportion
#>        1  0.3246   0.1362     0.3715     0.2715     0.1970 437.9     0.3079
#>        2 -0.7923  -0.6301    -0.9109    -0.7792    -0.9775 562.3     0.3954
#>        3  0.6925   0.6824     0.7982     0.7326     1.0747 421.8     0.2966
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```

In the staged solution, Profile 1 represents intermediate activity,
Profile 2 represents lower activity, and Profile 3 represents higher
activity. Their estimated `attendance` means are 0.197, -0.978, and
1.075, respectively. Although the numerical labels differ from those in
the earlier joint model, the same broad engagement patterns remain
recognizable. The estimated means also differ slightly because the
profiles were obtained under a different estimation procedure.

The parameters estimated at each stage can be examined through the
`"stages"` results table:

``` r

get_results(staged, "stages")
#>         stage group_classes            fixed log_likelihood parameters parameters_with_measurement
#> 1 measurement             1             <NA>          -8537         32                          32
#> 2  membership             3 means, variances          -8352          8                          38
#>   converged
#> 1      TRUE
#> 2      TRUE
```

The measurement stage estimates 32 parameters. The subsequent membership
stage estimates 8 additional parameters governing the student-class
composition. When the profile parameters retained from the measurement
stage are also counted, the membership-stage model contains 38
`parameters_with_measurement`. This latter quantity is useful when
making descriptive comparisons with a jointly estimated model because it
reflects both the newly estimated higher-level parameters and the
measurement parameters carried forward from the first stage.

The second-stage log likelihood is -8351.808, compared with -8348.429
for the joint model reported earlier. The lower log likelihood is
consistent with the additional restriction imposed by staged estimation:
the second stage cannot adjust the profile means and variances to
improve the overall fit because those parameters have already been
fixed. The two likelihoods should therefore be interpreted in light of
the different estimation constraints.

Overall, staged estimation provides a way to preserve a previously
established set of engagement profiles while estimating how those
profiles are distributed across student classes. Its main advantage in
this context is interpretive stability: the higher-level class structure
is estimated without allowing the course-level profile definitions to
change during the second stage.

## Interpretation and limitations

The membership regressions and grade contrasts indicate a consistent
relationship between prior academic performance and subsequent
engagement. Students in the lower-activity profile tend to have lower
previous grades than students in the more active profiles. The
distinction between the higher- and intermediate-activity profiles is
less clear, however, because these two profiles show greater
classification overlap and a correspondingly less certain difference in
previous grades.

The different estimation strategies should not be interpreted as
interchangeable versions of the same analysis. Three-step methods retain
the profile definitions established in the measurement model and then
examine their associations with external variables. In contrast,
one-step estimation allows the covariate to contribute directly to
estimation of the profiles and estimates the membership regression
conditional on student class. The choice between these approaches
therefore depends on the analytic objective: whether the priority is to
preserve an existing profile solution or to estimate the measurement and
covariate components jointly.

All of these results remain conditional on the adequacy of the
underlying measurement model. As discussed in the [evaluation
vignette](https://pak.dynasite.org/latents/articles/case-engagement-evaluation.md),
assumptions such as diagonal within-profile covariance structures, the
selected numbers of profiles and student classes, and the stability of
the likelihood solution require careful evaluation. Classification-error
corrections can account for uncertainty in assigning observations to
latent profiles, but they cannot compensate for a poorly specified
profile model.

Several inferential limitations also remain. Robust standard errors that
cluster enrolments within students account for dependence in the
observed data structure, but they do not convert an observational
association into a causal effect. In addition, both staged and
three-step procedures treat quantities estimated in earlier stages as
fixed. Their reported standard errors and confidence intervals therefore
do not propagate all of the uncertainty associated with estimating the
initial measurement model or classification-error matrix.
