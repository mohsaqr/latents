# Leisure activities of university students in daily life

Experience-sampling data on what university students did between prompts
when they were not studying. Prompts are nested within students, so the
data suit a two-level latent class analysis: each prompt belongs to an
activity profile, and students differ in how often their prompts fall in
each profile.

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

The students answered six prompts a day for fourteen days. At a prompt
where the student had not studied since the previous prompt, the
questionnaire asked which of eight leisure activities they had done. The
study asked these items only at non-study prompts, so the dataset
contains those prompts alone. Each row is one answered non-study prompt
with all eight activity items and all four affect ratings observed.

The dataset is a subset of the original study: 100 of its students,
drawn with a fixed seed from the students with at least 15 such prompts,
with every such prompt they answered. Each student contributes 15 to 51
prompts (median 24). The eight activity items are categorical
indicators; the four affect ratings are continuous-scored and can be
combined with them in a mixed measurement model. `worried` has a floor:
49.7% of its ratings are 1, so a profile of low-worry prompts can have
zero variance on it and reach the `min_variance` bound.

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
#> 1        1 time_with_friends       no  0.93676665  2.6956024
#> 2        2 time_with_friends       no  0.83220691  1.6013495
#> 3        1 time_with_friends      yes  0.06323335         NA
#> 4        2 time_with_friends      yes  0.16779309         NA
#> 5        1   on_social_media       no  0.33747581 -0.6745632
#> 6        2   on_social_media       no  0.92421483  2.5010419
#> 7        1   on_social_media      yes  0.66252419         NA
#> 8        2   on_social_media      yes  0.07578517         NA
#> 9        1    tv_video_games       no  0.63177676  0.5398467
#> 10       2    tv_video_games       no  0.86666716  1.8718065
#> 11       1    tv_video_games      yes  0.36822324         NA
#> 12       2    tv_video_games      yes  0.13333284         NA
#> 13       1    listened_music       no  0.64700412  0.6058963
#> 14       2    listened_music       no  0.94083907  2.7665107
#> 15       1    listened_music      yes  0.35299588         NA
#> 16       2    listened_music      yes  0.05916093         NA
#> 17       1            sports       no  0.94062608  2.7626906
#> 18       2            sports       no  0.91586401  2.3874334
#> 19       1            sports      yes  0.05937392         NA
#> 20       2            sports      yes  0.08413599         NA
#> 21       1           walking       no  0.88651172  2.0555948
#> 22       2           walking       no  0.83643759  1.6319573
#> 23       1           walking      yes  0.11348828         NA
#> 24       2           walking      yes  0.16356241         NA
#> 25       1           reading       no  0.87910684  1.9839993
#> 26       2           reading       no  0.89163129  2.1075133
#> 27       1           reading      yes  0.12089316         NA
#> 28       2           reading      yes  0.10836871         NA
#> 29       1     part_time_job       no  0.98376484  4.1042074
#> 30       2     part_time_job       no  0.95860865  3.1424110
#> 31       1     part_time_job      yes  0.01623516         NA
#> 32       2     part_time_job      yes  0.04139135         NA
get_results(lca, what = "profile_probabilities")
#>   group_class profile probability group_class_probability
#> 1           1       1   0.7093891               0.4830744
#> 2           1       2   0.2906109               0.4830744
#> 3           2       1   0.0762557               0.5169256
#> 4           2       2   0.9237443               0.5169256
# }
```
