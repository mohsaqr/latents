# Student engagement across a sequence of courses

A simulated dataset with the structure this package exists for, and all
of it in one place: enrolments nested in students, ordered within each
student, with a covariate that predicts which pattern an enrolment
follows.

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

  Integer course identifier, 1 to 32. Which course the enrolment was in.

- sequence:

  Integer, the position of this course in the student's own order, 1
  upwards. Pass it as `time`. Students who took fewer courses have
  shorter sequences, so the panel is ragged.

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

  Numeric, days the student was active in the course as on the
  `log1p(events)` scale, standardized within course. It is **generated
  from the four click measures** and then blurred, because a student is
  recorded present on a day precisely because they clicked something,
  and a day is a coarse unit – one tick however much happened inside it.
  The indicators are therefore locally dependent on purpose:
  `attendance` correlates 0.51 to 0.67 with the click measures overall
  and 0.16 to 0.26 within a single profile, which is more association
  than a shared profile alone would produce. See the note on local
  dependence below.

- previous_grade:

  Numeric, the grade the student earned in the course *before* this one,
  standardised across the cohort. At a student's first course an entry
  grade stands in its place. It is a predictor of engagement here, not
  an outcome of it: pass it as `profile_covariates` to
  `multilpa(profile_covariates = )`, or as the predictor to
  [`r3step()`](https://pak.dynasite.org/latents/reference/r3step.md).

- engagement:

  Factor, `"disengaged"` or `"engaged"`, the pattern the row was
  generated from. This is the truth behind the data, not something a
  study would observe; it is shipped so a fitted model can be checked
  against what produced it. It is deliberately not called `profile`:
  `get_results(fit, "assignments")` adds a column of that name, and a
  dataset already carrying it would make the verb refuse to join its own
  output.

- student_type:

  Factor, `"committed"` or `"wavering"`, the kind of student the row's
  student was generated as. The truth behind the group classes, shipped
  for the same reason.

## Source

Simulated by `data-raw/course-engagement.R`, seed 2026. The simulation
is calibrated to the shape of a real learning-analytics export – its
scales, its spreads, its correlations and how persistent engagement is
across a student's courses – but only those aggregate constants were
used. No row, identifier or value of any real student is present, and
the source is not distributed.

## Details

Each row is one student's enrolment in one course. A student takes up to
fifteen courses, and `sequence` orders them, so the same rows support a
cross-sectional two-level model (`id = "student"`), a latent transition
model (`time = "sequence"`), and a membership-covariate model
(`profile_covariates = "previous_grade"`).

Two kinds of student were generated. Both contain both patterns of
engagement; they differ in the *mix*, which is the two-level quantity
this package estimates. A committed student is engaged in about four
enrolments in five, a wavering student in about one in four. Pooled
across the sample the engaged share describes neither kind, which is the
point: recovering the split is what a two-level model does and a pooled
one cannot.

## Local independence does not hold here, on purpose

A latent profile model with a diagonal residual covariance assumes the
indicators are independent given the profile. `attendance` is built from
the click measures, so for these data that assumption is false, and the
package finds it: `get_results(fit, "residuals")` reports a largest
residual of 0.21 `p = 1.7e-15`, and `covariance_model = "full"` improves
the log likelihood from -7618.5 to -7494.7 for twenty more parameters,
which every information criterion prefers.

This is deliberate. A bundled dataset that satisfies local independence
exactly can only ever show those verbs agreeing that nothing is wrong;
here the diagnostic has a real violation to find, its cause is known,
and the fix can be seen to work. One consequence is worth knowing before
it surprises you: under a diagonal model
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
keeps preferring more classes, because extra classes absorb residual
association the model is not allowed to express directly.

The five activity indicators are **standardized within course**: each
value is that enrolment's distance from its own course's mean, in that
course's standard deviations. Every indicator therefore has mean zero
and standard deviation one within every course, and a `browse` of 1.2
reads as "1.2 standard deviations more active than the average enrolment
in that course". This is the scale a learning-analytics reader meets,
because activity volume is not comparable across courses: a course with
weekly quizzes generates more clicks than a seminar whatever its
students are like, so an analyst compares a student with their
coursemates rather than with the whole cohort.

In this simulation the courses do not in fact differ in activity volume,
so the standardization sets the scale rather than removing a course
effect. It is stated here rather than left implicit because a reader who
assumed otherwise would misread what the `course` column does.

The indicators are simulated directly on the `log1p(events)` scale, the
scale of a logged activity export; no raw counts are generated. They are
standardized on that scale, not as counts. Standardizing is affine
within a course, so it leaves the within-course skew exactly as it found
it, and raw learning-analytics counts are strongly right-skewed: on an
unlogged draft of these data the skew alone looked like an extra class,
exactly as a floor or ceiling effect does. The log removes the skew; the
standardization puts the courses on a common footing afterwards.

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
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.1703746               0.3252304
#> 2           1       2   0.8296254               0.3252304
#> 3           2       1   0.7858826               0.6747696
#> 4           2       2   0.2141174               0.6747696
# The generating truth ships alongside, so recovery can be checked directly.
get_results(fit, what = "assignments", data = course_engagement,
         truth = "engagement")
#>   assignment class      truth      value   n proportion
#> 1    profile     1 engagement disengaged  16 0.02735043
#> 2    profile     2 engagement disengaged 569 0.97264957
#> 3    profile     1 engagement    engaged 820 0.97968937
#> 4    profile     2 engagement    engaged  17 0.02031063
```
