# Simulated engagement states and trajectories over a study programme

`engagement_momentum` is a simulated dataset of student engagement
across the courses of a study programme. Each row represents one
student's enrolment in one course, with courses ordered within students.
In every course a student is in one of three engagement states (active,
average or disengaged), measured by five activity indicators
standardized within course; over the programme each student follows one
of three trajectories (engaged, fluctuating or mostly disengaged). Each
course also has a grade.

The design follows the study of engagement dynamics by López-Pernas and
Saqr (2024): engaged students show momentum, an active course tending to
be followed by another, and mostly disengaged students show inertia, a
disengaged course tending to be followed by another. The parameter
values are illustrative choices that reproduce these patterns; they are
not estimates from that study.

## Usage

``` r
engagement_momentum
```

## Format

A data frame with 1310 rows (106 students, 32 courses, one row per
student and course, 10 to 15 enrolments each, no missing values) and 12
columns:

- student:

  Integer student identifier, 1 to 106. The nesting unit; pass it as
  `id`.

- course:

  Integer course identifier, 1 to 32.

- sequence:

  Integer position of the course in the student's programme; pass it as
  `time`.

- browse:

  Numeric, course page views, standardized within course.

- lectures:

  Numeric, lecture views, standardized within course.

- forum_read:

  Numeric, forum posts read, standardized within course.

- forum_post:

  Numeric, forum posts written, standardized within course.

- active_days:

  Numeric, days with online activity, standardized within course.

- grade:

  Numeric course grade, standardized over the data.

- previous_grade:

  Numeric grade of the previous course (an entry grade at the first
  course), standardized over the data.

- state:

  Factor, the generating engagement state: `"active"`, `"average"` or
  `"disengaged"`.

- trajectory:

  Factor, the generating trajectory, constant within student:
  `"engaged"`, `"fluctuating"` or `"mostly_disengaged"`.

## Source

Simulated by `data-raw/engagement-momentum.R`, following the design of
López-Pernas, S., & Saqr, M. (2024). How the dynamics of engagement
explain the momentum of achievement and the inertia of disengagement: A
complex systems theory approach. *Computers in Human Behavior*, 153,
108126.
[doi:10.1016/j.chb.2023.108126](https://doi.org/10.1016/j.chb.2023.108126)

## Details

Students take 10 to 15 courses, ordered by `sequence`. Use
`id = "student"` for the nesting structure and `time = "sequence"` for
occasion order. The trajectories differ both in the mix of states they
visit, which
[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
recovers as group classes, and in their transitions between states,
which [`lta()`](https://pak.dynasite.org/latents/reference/lta.md)
recovers. The five indicators are drawn independently within a state, so
the default diagonal covariance structure is the generating one; unlike
`course_engagement`, the data contain no within-state dependence between
indicators. The grade of each course depends on its engagement state and
on the previous grade.

## Examples

``` r
descriptives(engagement_momentum,
             vars = c("browse", "lectures", "forum_read", "forum_post",
                      "active_days"),
             id = "student")
#>      variable    n n_missing          mean        sd   min  max n_distinct
#> 1      browse 1310         0 -8.396947e-05 0.9880070 -2.74 2.46        401
#> 2    lectures 1310         0 -5.343511e-05 0.9881515 -2.79 2.96        399
#> 3  forum_read 1310         0 -9.923664e-05 0.9880624 -2.46 2.80        386
#> 4  forum_post 1310         0  6.106870e-05 0.9879953 -2.54 3.06        395
#> 5 active_days 1310         0 -9.160305e-05 0.9880802 -2.62 2.83        392
#>   n_groups sd_within sd_between       icc
#> 1      106 0.8657215  0.5343503 0.2338112
#> 2      106 0.8847087  0.5087785 0.1998268
#> 3      106 0.8479255  0.5624311 0.2652751
#> 4      106 0.8781614  0.5179887 0.2114586
#> 5      106 0.8494755  0.5602438 0.2625956
fit <- multilpa(engagement_momentum,
                vars = c("browse", "lectures", "forum_read", "forum_post",
                         "active_days"),
                id = "student", n_profiles = 3, n_group_classes = 3,
                n_starts = 2, seed = 1)
get_results(fit, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1  0.85961489               0.1965273
#> 2           1       2  0.12216741               0.1965273
#> 3           1       3  0.01821770               0.1965273
#> 4           2       1  0.32321874               0.4461940
#> 5           2       2  0.34957303               0.4461940
#> 6           2       3  0.32720823               0.4461940
#> 7           3       1  0.06602054               0.3572788
#> 8           3       2  0.18266256               0.3572788
#> 9           3       3  0.75131689               0.3572788
# Compare the fitted classes with the generating states and trajectories.
get_results(fit, what = "assignments", data = engagement_momentum,
            truth = c("state", "trajectory"))
#>     assignment class      truth             value   n proportion
#> 1      profile     1      state            active   0 0.00000000
#> 2      profile     2      state            active  39 0.07103825
#> 3      profile     3      state            active 510 0.92896175
#> 4      profile     1      state           average  26 0.08049536
#> 5      profile     2      state           average 258 0.79876161
#> 6      profile     3      state           average  39 0.12074303
#> 7      profile     1      state        disengaged 415 0.94748858
#> 8      profile     2      state        disengaged  23 0.05251142
#> 9      profile     3      state        disengaged   0 0.00000000
#> 10 group_class     1 trajectory           engaged   0 0.00000000
#> 11 group_class     2 trajectory           engaged   6 0.13953488
#> 12 group_class     3 trajectory           engaged  37 0.86046512
#> 13 group_class     1 trajectory       fluctuating   1 0.02439024
#> 14 group_class     2 trajectory       fluctuating  39 0.95121951
#> 15 group_class     3 trajectory       fluctuating   1 0.02439024
#> 16 group_class     1 trajectory mostly_disengaged  20 0.90909091
#> 17 group_class     2 trajectory mostly_disengaged   2 0.09090909
#> 18 group_class     3 trajectory mostly_disengaged   0 0.00000000
```
