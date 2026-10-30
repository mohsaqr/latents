# Student engagement across a sequence of courses

`course_engagement` is a simulated dataset for illustrating student
engagement across successive courses. Each row represents one student's
enrolment in one course, with courses ordered within students. The data
include five activity indicators, standardised within course, and a
simulated previous-course grade.

The data were generated from two engagement profiles and two
student-level classes that differ in the prevalence of those profiles.
Engagement can change across a student's course sequence. The generating
profile and student-class labels are retained, allowing fitted
classifications to be compared with the structure used to generate the
data.

The dataset supports examples of multilevel profile analysis, latent
transition analysis and analyses relating previous grades to profile
membership.

## Usage

``` r
course_engagement
```

## Format

A data frame with 1422 rows (106 students, 32 courses, up to 15
enrolments each) and 11 columns:

- student:

  Integer student identifier, 1 to 106. The nesting unit; pass it as
  `id`.

- course:

  Integer course identifier, 1 to 32.

- sequence:

  Integer position of the course in the student's own order, 1 to 15.
  Pass it as `time`. Students who took fewer courses have shorter
  sequences, so the panel is unbalanced.

- browse:

  Numeric, simulated browsing activity, standardized within course.

- lectures:

  Numeric, simulated lecture-viewing activity, standardized within
  course.

- forum_read:

  Numeric, simulated forum-reading activity, standardized within course.

- forum_post:

  Numeric, simulated forum-posting activity, standardized within course.

- attendance:

  Numeric, simulated attendance activity, standardized within course.
  Generated from the four other activity indicators with added noise,
  introducing within-profile dependence; see *Local dependence*.

- previous_grade:

  Numeric, simulated grade in the previous course, standardized across
  the cohort; an entry grade for the first course. Can be used as a
  membership predictor through the `profile_covariates` argument of
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or in
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md).

- engagement:

  Factor, `"disengaged"` or `"engaged"`: the pattern used to generate
  the observation; retained for assessing class recovery.

- student_type:

  Factor, `"committed"` or `"wavering"`: the student's generating class;
  retained for assessing group-class recovery.

## Source

Simulated by `data-raw/course-engagement.R` in the source repository,
seed 2026. The simulation is calibrated to aggregate constants of a real
learning-analytics export (scales, spreads, correlations and the
persistence of engagement across courses). No row, identifier or value
of any real student is included, and the source data are not
distributed.

## Details

Students take 10 to 15 courses, ordered by `sequence`. Use
`id = "student"` for the nesting structure, `time = "sequence"` for
occasion order, and `profile_covariates = "previous_grade"` for a
membership-covariate model.

In the generating model, students in the `"committed"` class have a
stationary probability of an engaged enrolment of approximately 0.8,
compared with approximately 0.25 in the `"wavering"` class.

The activity indicators were simulated as continuous values intended to
represent log-transformed activity measures, then standardized within
course and rounded. Raw event counts were not generated. Each indicator
has approximately mean zero and standard deviation one within course;
for example, `browse = 1.2` represents a value 1.2 standard deviations
above the course mean. The simulation does not include between-course
differences in activity volume.

## Local dependence

`attendance` is generated from the four click measures, introducing
within-profile dependence among the indicators. A model with diagonal
residual covariance therefore does not reproduce the generating
covariance structure. Use `get_results(fit, "residuals")` to examine
residual associations and `covariance_model = "full"` to estimate
within-profile covariances. Class enumeration under a diagonal model may
represent residual association by additional classes. The two generating
profiles therefore do not imply that every fitted specification will
favour two profiles.

## See also

[`vignette("lpa", package = "latents")`](https://pak.dynasite.org/latents/articles/lpa.md)
for the analysis these data are used in,
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for the two-level model,
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) for the
sequence, and
[`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md) for
membership-covariate analysis.

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
#> 2        1   lectures  0.4274031 0.8446923          0.9190714
#> 3        1 forum_read  0.6198249 0.4814743          0.6938835
#> 4        1 forum_post  0.5304888 0.6887895          0.8299335
#> 5        1 attendance  0.6490498 0.3841263          0.6197793
#> 6        2     browse -0.7655998 0.5284684          0.7269583
#> 7        2   lectures -0.6064852 0.5383733          0.7337392
#> 8        2 forum_read -0.8792137 0.3631130          0.6025886
#> 9        2 forum_post -0.7526883 0.4211624          0.6489703
#> 10       2 attendance -0.9206599 0.3736293          0.6112523
#>    mean_standard_error variance_standard_error
#> 1           0.02710160              0.02948401
#> 2           0.03248649              0.04210054
#> 3           0.02473030              0.02439626
#> 4           0.02953027              0.03481645
#> 5           0.02226482              0.01955858
#> 6           0.03121144              0.03221700
#> 7           0.03093059              0.03209025
#> 8           0.02621510              0.02279295
#> 9           0.02768122              0.02607600
#> 10          0.02642702              0.02278996
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.1704508                0.325297
#> 2           1       2   0.8295492                0.325297
#> 3           2       1   0.7859317                0.674703
#> 4           2       2   0.2140683                0.674703
# Compare fitted assignments with the generating engagement profiles.
get_results(fit, what = "assignments", data = course_engagement,
         truth = "engagement")
#>   assignment class      truth      value   n proportion
#> 1    profile     1 engagement disengaged  16 0.02735043
#> 2    profile     2 engagement disengaged 569 0.97264957
#> 3    profile     1 engagement    engaged 820 0.97968937
#> 4    profile     2 engagement    engaged  17 0.02031063
```
