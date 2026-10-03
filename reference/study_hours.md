# Study hours and quiz scores under two study strategies

A simulated dataset for
[`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md),
built so that every nesting of the mixture regression has a model that
generated some of its columns, and the truth is kept for checking
recovery.

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

  Factor, `"deep"` or `"surface"`: the week's generating strategy. Not
  for fitting; kept to check recovery.

- student_type:

  Factor, `"steady"` or `"erratic"`: the student's generating kind. Not
  for fitting; kept to check recovery.

## Source

Simulated by `data-raw/study-hours.R` in the source repository.

## Details

150 students each report six study weeks. In every week a student
studies with one of two strategies, and the strategy decides how hours
become a quiz `score`: a deep week gains about 4.5 points per hour from
a base of 35, a surface week under one point per hour from a base of 55.
Students are of two kinds that differ in their *mix* of strategies –
steady students study deeply in about 85% of weeks, erratic ones in
about 25% – and higher `motivation` makes a student more likely to be
steady. More `sleep` makes a deep week more likely.

So `score ~ hours` is a two-class mixture at the level of the week
(`class_level = "observation"`), whose class shares differ by student
(the two-level model, `id = "student"`, `n_group_classes = 2`). `passed`
depends on the student's kind, not the week's strategy, so it is a
mixture at the level of the student (`class_level = "group"`).
`questions` is a count whose Poisson regression depends on the week's
strategy.

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
#> Class 1  Intercept     34.35  [32.94, 35.76]  <.001
#> Class 1  hours          4.57  [ 4.36,  4.77]  <.001
#> Class 2  Intercept     54.90  [52.75, 57.04]  <.001
#> Class 2  hours          0.81  [ 0.48,  1.14]  <.001
```
