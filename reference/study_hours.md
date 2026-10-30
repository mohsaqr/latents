# Study hours and quiz scores under two study strategies

`study_hours` is a simulated dataset illustrating how the relationship
between study time and quiz performance can differ across latent
classes. Each row represents one study week for a student, with six
weeks recorded for each of 150 students.

The data were generated from two week-level regression classes, labelled
"deep" and "surface", with different intercepts and slopes relating
study hours to quiz scores. Two student-level classes differ in the
prevalence of these weekly classes, allowing class membership to vary
across weeks within students.

The generating labels are retained for assessing class recovery.
Additional binary and count outcomes provide examples of mixture
regression with other response distributions.

## Usage

``` r
study_hours
```

## Format

A data frame with 900 rows (150 students by 6 weeks) and 10 columns:

- student:

  Integer student identifier, 1 to 150.

- week:

  Integer study week, 1 to 6.

- hours:

  Hours studied that week, 0 to 12.

- sleep:

  Hours slept the night before the quiz.

- motivation:

  Standardized motivation, constant within student.

- score:

  Quiz score.

- passed:

  Integer 0/1: passed the weekly check.

- questions:

  Count of questions the student asked that week.

- strategy:

  Factor, `"deep"` or `"surface"`: the week's generating strategy;
  retained for assessing class recovery.

- student_type:

  Factor, `"steady"` or `"erratic"`: the student's generating class;
  retained for assessing class recovery.

## Source

Simulated by `data-raw/study-hours.R` in the source repository.

## Details

The generating model contains two week-level study strategies with
different regressions of quiz `score` on `hours`. The `"deep"` strategy
has an intercept of approximately 35 and a slope of 4.5 points per hour;
the `"surface"` strategy has an intercept of approximately 55 and a
slope below one point per hour. Students in the `"steady"` class use the
deep strategy in approximately 85% of weeks, compared with 25% for the
`"erratic"` class. In the simulation, `motivation` predicts
student-level membership and `sleep` predicts week-level membership.

`score ~ hours` is therefore a two-class mixture at the week level
(`class_level = "observation"`) with class shares that differ by student
(`id = "student"`, `n_group_classes = 2`). The binary outcome `passed`
depends on the student's generating class and illustrates a mixture at
the student level (`class_level = "group"`). `questions` is a count
whose Poisson regression depends on the week's strategy.

## See also

[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
[`enumerate_regressions()`](https://pak.dynasite.org/latents/reference/enumerate_regressions.md).

## Examples

``` r
fit <- mixture_regression(score ~ hours, data = study_hours, n_classes = 2,
                          n_starts = 3, seed = 1)
get_results(fit, "coefficients")
#> Regression coefficients (95% CI)
#> 
#> Class    Term       Estimate          95% CI      p
#> -------  ---------  --------  --------------  -----
#> Class 1  Intercept     34.36  [32.95, 35.77]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.76, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
```
