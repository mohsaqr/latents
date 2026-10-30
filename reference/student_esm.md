# Leisure activities of university students in daily life

`student_esm` contains experience-sampling reports of leisure activities
and affect from 100 university students. Each row represents a survey
response collected when the student reported no studying since the
previous prompt. Eight binary items describe activities undertaken
during that interval, and four ratings describe affect.

Repeated reports are nested within students, providing an example for
multilevel latent class analysis of activity patterns and differences in
their prevalence between students. The dataset retains complete reports
from students with at least 15 eligible responses. The affect ratings
can also be included as indicators in a mixed measurement model.

## Usage

``` r
student_esm
```

## Format

A data frame with 2582 rows (100 students) and 15 columns:

- student:

  Integer student identifier, 1 to 100. The nesting unit; pass it as
  `id`.

- day:

  Integer study day, 0 to 13.

- beep:

  Integer prompt of the day, 1 to 5 among the prompts kept.

- time_with_friends:

  Factor, `"no"` or `"yes"`: spent time with friends since the previous
  prompt.

- on_social_media:

  Factor, `"no"` or `"yes"`: used social media.

- tv_video_games:

  Factor, `"no"` or `"yes"`: watched TV or played video games.

- listened_music:

  Factor, `"no"` or `"yes"`: listened to music.

- sports:

  Factor, `"no"` or `"yes"`: did sports.

- walking:

  Factor, `"no"` or `"yes"`: went for a walk.

- reading:

  Factor, `"no"` or `"yes"`: read.

- part_time_job:

  Factor, `"no"` or `"yes"`: worked in a part-time job.

- happy:

  Integer rating of feeling happy, 1 to 7.

- relaxed:

  Integer rating of feeling relaxed, 1 to 7.

- worried:

  Integer rating of feeling worried, 1 to 7.

- exhausted:

  Integer rating of feeling exhausted, 1 to 7.

## Source

openESM dataset 0062, deposited by Neubauer and Schmiedek at Zenodo,
[doi:10.5281/zenodo.17347974](https://doi.org/10.5281/zenodo.17347974) ,
under the Creative Commons Attribution 4.0 licence (CC-BY 4.0), and
retrieved through the openESM database (<https://openesmdata.org>). The
subset is produced by `data-raw/student-esm.R` in the source repository;
identifiers are renumbered and the values are otherwise unchanged.

## Details

In the original study students answered six prompts a day for fourteen
days; the activity items were asked only at prompts where the student
had not studied since the previous prompt. The dataset contains 100
students, drawn with a fixed seed from those with at least 15 complete
non-study reports, retaining every such report: 15 to 51 per student
(median 24). The eight activity items are categorical indicators; the
four affect ratings can be added as continuous indicators in a mixed
measurement model. `worried` has a pronounced floor effect (49.7% of
ratings are 1). A Gaussian component representing low-worry prompts may
therefore reach the `min_variance` bound.

## References

Neubauer, A. B., & Schmiedek, F. (2024). Approaching academic adjustment
on multiple time scales. *Zeitschrift für Erziehungswissenschaft*, 27,
147–168.
[doi:10.1007/s11618-023-01182-8](https://doi.org/10.1007/s11618-023-01182-8)

## See also

[`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
with `categorical =` for the two-level latent class model.

## Examples

``` r
head(student_esm)
#>   student day beep time_with_friends on_social_media tv_video_games
#> 1       1   0    3               yes             yes            yes
#> 2       1   0    5                no              no            yes
#> 3       1   1    1                no              no             no
#> 4       1   1    3                no              no            yes
#> 5       1   1    4                no              no            yes
#> 6       1   2    3                no              no             no
#>   listened_music sports walking reading part_time_job happy relaxed worried
#> 1             no    yes      no     yes            no     5       6       5
#> 2             no     no      no      no            no     5       6       2
#> 3             no     no      no      no            no     3       3       2
#> 4             no     no      no      no            no     6       5       2
#> 5             no     no      no      no            no     4       5       1
#> 6             no     no      no      no            no     5       4       5
#>   exhausted
#> 1         3
#> 2         5
#> 3         7
#> 4         4
#> 5         6
#> 6         7
# \donttest{
activities <- c("time_with_friends", "on_social_media", "tv_video_games",
                "listened_music", "sports", "walking", "reading",
                "part_time_job")
first_week <- subset(student_esm, day <= 6)
lca <- multilpa(first_week, vars = activities, id = "student",
                n_profiles = 2, n_group_classes = 2,
                categorical = activities, n_starts = 3, seed = 1)
get_results(lca, what = "responses")
#>    profile         indicator category probability  threshold
#> 1        1 time_with_friends       no  0.93678831  2.6959681
#> 2        2 time_with_friends       no  0.83243250  1.6029659
#> 3        1 time_with_friends      yes  0.06321169         NA
#> 4        2 time_with_friends      yes  0.16756750         NA
#> 5        1   on_social_media       no  0.33643583 -0.6792181
#> 6        2   on_social_media       no  0.92352634  2.4912530
#> 7        1   on_social_media      yes  0.66356417         NA
#> 8        2   on_social_media      yes  0.07647366         NA
#> 9        1    tv_video_games       no  0.63130098  0.5378021
#> 10       2    tv_video_games       no  0.86642891  1.8697462
#> 11       1    tv_video_games      yes  0.36869902         NA
#> 12       2    tv_video_games      yes  0.13357109         NA
#> 13       1    listened_music       no  0.64636806  0.6031124
#> 14       2    listened_music       no  0.94056673  2.7616284
#> 15       1    listened_music      yes  0.35363194         NA
#> 16       2    listened_music      yes  0.05943327         NA
#> 17       1            sports       no  0.94060089  2.7622397
#> 18       2            sports       no  0.91593649  2.3883744
#> 19       1            sports      yes  0.05939911         NA
#> 20       2            sports      yes  0.08406351         NA
#> 21       1           walking       no  0.88653599  2.0558360
#> 22       2           walking       no  0.83653689  1.6326833
#> 23       1           walking      yes  0.11346401         NA
#> 24       2           walking      yes  0.16346311         NA
#> 25       1           reading       no  0.87896299  1.9826464
#> 26       2           reading       no  0.89169308  2.1081529
#> 27       1           reading      yes  0.12103701         NA
#> 28       2           reading      yes  0.10830692         NA
#> 29       1     part_time_job       no  0.98375148  4.1033718
#> 30       2     part_time_job       no  0.95867460  3.1440743
#> 31       1     part_time_job      yes  0.01624852         NA
#> 32       2     part_time_job      yes  0.04132540         NA
#>    probability_standard_error
#> 1                 0.012814476
#> 2                 0.015246632
#> 3                 0.012814476
#> 4                 0.015246632
#> 5                 0.045916880
#> 6                 0.026634764
#> 7                 0.045916880
#> 8                 0.026634764
#> 9                 0.028916367
#> 10                0.016049952
#> 11                0.028916367
#> 12                0.016049952
#> 13                0.031673364
#> 14                0.014195450
#> 15                0.031673364
#> 16                0.014195450
#> 17                0.012077514
#> 18                0.010533067
#> 19                0.012077514
#> 20                0.010533067
#> 21                0.016567397
#> 22                0.014229716
#> 23                0.016567397
#> 24                0.014229716
#> 25                0.016848516
#> 26                0.011899028
#> 27                0.016848516
#> 28                0.011899028
#> 29                0.006254343
#> 30                0.007299599
#> 31                0.006254343
#> 32                0.007299599
get_results(lca, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1  0.07511869               0.5171722
#> 2           1       2  0.92488131               0.5171722
#> 3           2       1  0.70803548               0.4828278
#> 4           2       2  0.29196452               0.4828278
# }
```
