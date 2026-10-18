# Student engagement across a sequence of courses

A simulated dataset of course enrolments nested in students. Each row is
one student's enrolment in one course. Students take 10 to 15 courses,
ordered by `sequence`, so the data support a cross-sectional two-level
model (`id = "student"`), a latent transition model
(`time = "sequence"`) and a membership-covariate model
(`profile_covariates = "previous_grade"`).

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

  Numeric, course page views on the `log1p(events)` scale, standardized
  within course.

- lectures:

  Numeric, lecture videos viewed on the `log1p(events)` scale,
  standardized within course.

- forum_read:

  Numeric, forum posts read on the `log1p(events)` scale, standardized
  within course.

- forum_post:

  Numeric, forum posts written on the `log1p(events)` scale,
  standardized within course.

- attendance:

  Numeric, days the student was active in the course, on the
  `log1p(events)` scale, standardized within course. Generated from the
  four click measures with added noise, so it correlates 0.51 to 0.67
  with them overall and 0.10 to 0.23 among the enrolments assigned to
  one profile of the example fit. See *Local dependence*.

- previous_grade:

  Numeric, the grade the student earned in the previous course,
  standardised across the cohort; an entry grade at a student's first
  course. A predictor of engagement, for
  `multilpa(profile_covariates = )` or
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md).

- engagement:

  Factor, `"disengaged"` or `"engaged"`: the pattern the row was
  generated from. Not for fitting; kept to check recovery.

- student_type:

  Factor, `"committed"` or `"wavering"`: the student's generating kind,
  corresponding to the group classes. Not for fitting; kept to check
  recovery.

## Source

Simulated by `data-raw/course-engagement.R`, seed 2026. The simulation
is calibrated to aggregate constants of a real learning-analytics export
(scales, spreads, correlations and the persistence of engagement across
courses). No row, identifier or value of any real student is included,
and the source data are not distributed.

## Details

Enrolments follow one of two engagement patterns, and students are of
two kinds that differ in the mix of patterns: a committed student is
engaged in about four enrolments in five, a wavering student in about
one in four. The pooled engaged share describes neither kind; a
two-level model recovers the two kinds.

The five activity indicators were simulated on the `log1p(events)` scale
of a logged activity export (no raw counts were generated) and then
standardised within course: each value is the enrolment's distance from
its course's mean in that course's standard deviations, so every
indicator has mean zero and standard deviation one within every course.
A `browse` of 1.2 is 1.2 standard deviations above the average enrolment
in that course. The simulated courses do not differ in activity volume,
so the standardisation sets the scale but removes no course effect. The
log transformation removes the right skew of raw activity counts, which
standardisation within course would leave unchanged.

## Local dependence

`attendance` is generated from the four click measures, so the
indicators are not independent given the profile, as a latent profile
model with a diagonal residual covariance assumes. For the two-profile,
two-group-class fit of the example, `get_results(fit, "residuals")`
reports a largest residual of 0.23 (`lectures` with `attendance` within
the first profile, `p = 6.4e-12`), and `covariance_model = "full"`
improves the log likelihood from -8439.8 to -8313.1 for twenty more
parameters, preferred by every information criterion except the
individual-level AWE. Under a diagonal model,
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
tends to favour more classes, which absorb the residual association.

## See also

[`vignette("lpa", package = "latents")`](https://pak.dynasite.org/latents/articles/lpa.md)
for the analysis these data are used in,
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
for the two-level model,
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md) for the
sequence, and `multilpa(profile_covariates = )` for `previous_grade`.

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
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.7859317               0.6747028
#> 2           1       2   0.2140683               0.6747028
#> 3           2       1   0.1704509               0.3252972
#> 4           2       2   0.8295491               0.3252972
# The generating truth ships alongside, so recovery can be checked directly.
get_results(fit, what = "assignments", data = course_engagement,
         truth = "engagement")
#>   assignment class      truth      value   n proportion
#> 1    profile     1 engagement disengaged  16 0.02735043
#> 2    profile     2 engagement disengaged 569 0.97264957
#> 3    profile     1 engagement    engaged 820 0.97968937
#> 4    profile     2 engagement    engaged  17 0.02031063
```
