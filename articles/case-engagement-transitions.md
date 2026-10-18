# Case study: Latent transition analysis of course engagement

Latent transition analysis (LTA) estimates the probability of moving
from one latent profile to another between consecutive occasions
(Collins & Lanza, 2010). With group classes, each class has its own
starting distribution and its own transition matrix, so the classes
separate students who differ in how their engagement changes from
students who differ only in where it starts.

## Data and model

`course_engagement` holds 1,422 course enrolments of 106 students,
ordered within each student by `sequence`. The five indicators are
simulated as continuous values representing a log-count scale and
standardized within course; no raw counts or
[`log1p()`](https://rdrr.io/r/base/Log.html) transformation are
generated. The data are synthetic and revisit the design of Saqr and
López-Pernas (2021), who asked how students’ online engagement develops
across a program. That study identified course-level engagement states
and longer-term trajectories with latent class analysis and hidden
Markov models. Here, three profiles and two student classes are working
choices for illustrating the transition model; they do not establish the
number of profiles or reproduce the published estimates.

For student $`i`$ at course position $`t`$, $`C_{it}`$ is the latent
profile and $`G_i`$ the student’s class. The model estimates

- the initial probabilities $`\delta_{k\mid m}=P(C_{i1}=k\mid G_i=m)`$,
- the transition probabilities
  $`a_{kl\mid m}=P(C_{i,t+1}=l\mid C_{it}=k,G_i=m)`$, which sum to one
  across destinations $`l`$,
- and a Gaussian mean and variance for each profile and indicator.

The transitions are first order and homogeneous: the next profile
depends only on the current profile and the student’s class, and the
probabilities are the same at every position. Profile means and
variances are shared across positions and classes, so a change of
profile is a change of state. Indicators are independent within a
profile under the default diagonal covariance. The fit can fail to reach
the global maximum from a single start, and a transition that is rare
within a class is estimated from little information; the expected count
of each transition shows how much.

``` r

library(latents)
vars <- c("browse", "lectures", "forum_read", "forum_post", "attendance")

moves <- lta(course_engagement, vars, "student", n_profiles = 3,
             n_group_classes = 2, time = "sequence", n_starts = 10, seed = 1)
get_results(moves, "starts")
#>    start log_likelihood converged iterations error
#> 1      1          -8243      TRUE        133  <NA>
#> 2      2          -8243      TRUE        158  <NA>
#> 3      3          -8243      TRUE        132  <NA>
#> 4      4          -8243      TRUE        188  <NA>
#> 5      5          -8243      TRUE        116  <NA>
#> 6      6          -8243      TRUE        326  <NA>
#> 7      7          -8243      TRUE        143  <NA>
#> 8      8          -8243      TRUE        165  <NA>
#> 9      9          -8243      TRUE        135  <NA>
#> 10    10          -8243      TRUE        160  <NA>
```

All ten starts converge to the same log likelihood, so the solution is
replicated within this seed’s starts.

## The profiles

[`plot_profiles()`](https://pak.dynasite.org/latents/reference/latents-plots.md)
draws the profile means across the indicators, and each line’s label
gives the profile’s share of the enrolments.

``` r

plot_profiles(moves)
```

![](case-engagement-transitions_files/figure-html/profiles-1.png)

``` r

get_results(moves)
#> No standard errors in this table: Wald inference is unavailable when a probability sits at its lower bound: its logit is at minus infinity, where the likelihood has no curvature to invert.
#>    profile  indicator    mean variance standard_deviation mean_standard_error
#> 1        1     browse -0.7760   0.5149             0.7175                  NA
#> 2        1   lectures -0.6117   0.5417             0.7360                  NA
#> 3        1 forum_read -0.8881   0.3540             0.5950                  NA
#> 4        1 forum_post -0.7607   0.4149             0.6441                  NA
#> 5        1 attendance -0.9293   0.3727             0.6105                  NA
#> 6        2     browse  0.7058   0.5292             0.7274                  NA
#> 7        2   lectures  0.7115   0.7967             0.8926                  NA
#> 8        2 forum_read  0.8070   0.4653             0.6821                  NA
#> 9        2 forum_post  0.7469   0.5767             0.7594                  NA
#> 10       2 attendance  1.1240   0.2067             0.4546                  NA
#> 11       3     browse  0.3959   0.6100             0.7810                  NA
#> 12       3   lectures  0.1892   0.7592             0.8713                  NA
#> 13       3 forum_read  0.4541   0.4578             0.6766                  NA
#> 14       3 forum_post  0.3445   0.7146             0.8453                  NA
#> 15       3 attendance  0.2538   0.2044             0.4521                  NA
#>    variance_standard_error
#> 1                       NA
#> 2                       NA
#> 3                       NA
#> 4                       NA
#> 5                       NA
#> 6                       NA
#> 7                       NA
#> 8                       NA
#> 9                       NA
#> 10                      NA
#> 11                      NA
#> 12                      NA
#> 13                      NA
#> 14                      NA
#> 15                      NA
```

Profile 1 is the lowest on every indicator, profile 2 the highest, and
profile 3 lies between them; their attendance means are -0.929, 1.124
and 0.254. The profiles are therefore lower, higher and intermediate
activity. The heatmap shows the same means as standard deviations from
each indicator’s grand mean.

``` r

plot_profiles(moves, what = "heatmap")
```

![](case-engagement-transitions_files/figure-html/heatmap-1.png)

The contrast between profiles is present on all five indicators, so the
profiles differ in overall activity more than in its composition.

## Where students start

In the `"initial"` table, `probability` is $`\delta_{k\mid m}`$, the
chance of starting in each profile. `prevalence` is the expected share
of all of a class’s enrolments in that profile, across the whole
history, and `group_class_probability` is the share of students in the
class.

``` r

get_results(moves, "initial")
#>   group_class profile probability prevalence group_class_probability
#> 1           1       1   2.037e-01    0.22973                  0.6996
#> 2           1       2   4.437e-01    0.33313                  0.6996
#> 3           1       3   3.526e-01    0.43713                  0.6996
#> 4           2       1   7.623e-01    0.81929                  0.3004
#> 5           2       2   2.377e-01    0.09657                  0.3004
#> 6           2       3   6.814e-07    0.08414                  0.3004
```

Student class 2, with 30.0% of students, starts in the lower-activity
profile with probability 0.762 and in the higher-activity profile with
probability 0.238; its probability of starting in the intermediate
profile is effectively zero. Class 1, with 70.0% of students, starts
across all three profiles, with probabilities 0.204, 0.444 and 0.353.
Starting and overall shares differ: class 2 starts in the lower profile
with probability 0.762, but 0.819 of its enrolments are expected there
over the whole history, while class 1 starts in the higher profile with
probability 0.444 and spends 0.333 of its enrolments there.

## How engagement changes

The `"transitions"` table gives $`a_{kl\mid m}`$ for every class and
pair of profiles. `stable` marks staying in the same profile.
`expected_count` is the posterior expected number of such transitions in
the data, so it can be fractional, and it shows how much information
supports each estimate. `estimated = FALSE` would mark a row the data
cannot inform, because the profile is never occupied before a student’s
last course; such a row holds uniform placeholder values. Every row is
estimated in this fit.

``` r

get_results(moves, "transitions")
#>    group_class from to probability expected_count stable estimated group_class_probability
#> 1            1    1  1     0.61456        127.937   TRUE      TRUE                  0.6996
#> 2            1    1  2     0.17562         36.559  FALSE      TRUE                  0.6996
#> 3            1    1  3     0.20982         43.680  FALSE      TRUE                  0.6996
#> 4            1    2  1     0.12770         38.858  FALSE      TRUE                  0.6996
#> 5            1    2  2     0.37733        114.816   TRUE      TRUE                  0.6996
#> 6            1    2  3     0.49496        150.613  FALSE      TRUE                  0.6996
#> 7            1    3  1     0.11407         46.483  FALSE      TRUE                  0.6996
#> 8            1    3  2     0.36049        146.893  FALSE      TRUE                  0.6996
#> 9            1    3  3     0.52544        214.119   TRUE      TRUE                  0.6996
#> 10           2    1  1     0.93856        303.231   TRUE      TRUE                  0.3004
#> 11           2    1  2     0.03126         10.101  FALSE      TRUE                  0.3004
#> 12           2    1  3     0.03017          9.748  FALSE      TRUE                  0.3004
#> 13           2    2  1     0.29609         11.618  FALSE      TRUE                  0.3004
#> 14           2    2  2     0.48188         18.909   TRUE      TRUE                  0.3004
#> 15           2    2  3     0.22203          8.712  FALSE      TRUE                  0.3004
#> 16           2    3  1     0.33924         11.440  FALSE      TRUE                  0.3004
#> 17           2    3  2     0.14060          4.741  FALSE      TRUE                  0.3004
#> 18           2    3  3     0.52016         17.541   TRUE      TRUE                  0.3004
```

The transition plot draws each class’s matrix as a heatmap, with the
origin profile in rows and the destination in columns.

``` r

plot_transitions(moves)
```

![](case-engagement-transitions_files/figure-html/transition-plot-1.png)

In class 2, lower activity persists with probability 0.939, estimated
from 303 expected transitions. Each of its other transitions rests on
fewer than 20 expected transitions, because class-2 students rarely
leave the lower profile. Class 1 moves between profiles. Its students
remain in the lower profile with probability 0.615, move from higher to
intermediate activity with probability 0.495, and from intermediate to
higher with probability 0.360. Class 2 is therefore a persistently
low-engagement class, and class 1 a class whose engagement changes from
course to course.

## Sparse and boundary estimates

Some estimates rest on very little information. Class 2’s probability of
starting in the intermediate profile is at its lower bound, because no
class-2 student is expected to start there. Class 2’s transitions out of
the higher and intermediate profiles rest on between 4.7 and 18.9
expected transitions each, because its students rarely leave the lower
profile. A probability near zero shows that the data carry little
information about that movement, not that the movement is impossible.
Estimates like these can also mean that the chosen number of classes
divides the available transitions too finely.

## Assigned sequences

The sequence plot shows each student’s modal profile at each course
position, one row per student. Modal assignment discards classification
uncertainty.

``` r

plot_sequences(moves)
```

![](case-engagement-transitions_files/figure-html/sequences-1.png)

``` r

get_results(moves, "sequence_summary")
#>   group_class groups observations mean_length median_length shortest longest complete gaps
#> 1           1     75         1005       13.40            14       10      15       16    0
#> 2           2     31          417       13.45            13       10      15        8    0
```

The summary gives the number of students and courses in each assigned
class. Sequence lengths differ between students, so the rows of the plot
end at different positions.

## Classification uncertainty

``` r

get_results(moves, "entropy")
#>         level n_classes n_units entropy_sum relative_entropy
#> 1 individuals         3    1422      309.45           0.8019
#> 2      groups         2     106       16.31           0.7781
```

Relative entropy is 0.802 for enrolments and 0.778 for students. The
entropy plot shows, within each profile, how the per-case entropy
contributions are distributed, and the posterior plot shows the
posterior probability of each case’s assigned profile.

``` r

plot_diagnostics(moves, what = "entropy")
```

![](case-engagement-transitions_files/figure-html/entropy-plot-1.png)

``` r

plot_diagnostics(moves, what = "posteriors")
```

![](case-engagement-transitions_files/figure-html/posteriors-plot-1.png)

Both levels retain uncertainty, so an analysis that replaces the
posteriors by modal assignments treats some enrolments and students as
more certain than the model estimates them to be.

## Transition networks

[`get_tna()`](https://pak.dynasite.org/latents/reference/get_tna.md)
returns the marginal transition matrix as a model of the tna package
(Tikka et al., 2025). It pools the posterior expected transition counts
across classes and normalizes them, so it uses the estimated
transitions, not transitions counted from modal assignments.
[`get_group_tna()`](https://pak.dynasite.org/latents/reference/get_group_tna.md)
returns one network per class.

``` r

network <- get_tna(moves)
plot(network)
```

![](case-engagement-transitions_files/figure-html/network-1.png)

``` r

per_class <- get_group_tna(moves)
plot(per_class)
```

![](case-engagement-transitions_files/figure-html/group-network-1.png)![](case-engagement-transitions_files/figure-html/group-network-2.png)

The class networks show the two patterns of the transition table: a
dominant self-loop on the lower profile in class 2, and movement among
all three profiles in class 1.

## Limitations

The model assumes constant transition probabilities across course
positions, so change specific to one stage of the program is averaged
over positions. Transitions are between successive observed courses,
which can differ in duration, content and assessment. Class differences
describe patterns in engagement histories and do not identify their
causes. Wald standard errors are withheld for this fit, because class
2’s probability of starting in the intermediate profile sits at its
lower bound (the note above the profile table), and the numbers of
profiles and classes were not compared here;
[`enumerate_classes()`](https://pak.dynasite.org/latents/reference/enumerate_classes.md)
with `time` compares transition models over a grid of profile and class
counts. The number of profiles and classes here is therefore a working
choice. The diagonal covariance ignores the association between
`attendance` and the click measures, which the evaluation vignette shows
is present in these data; unmodelled dependence can be absorbed by extra
profiles and change the apparent transitions. Transitions supported by
few expected counts, such as those out of the higher and intermediate
profiles in class 2, should not be interpreted.

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent class and latent
transition analysis*. Wiley.

Saqr, M., & López-Pernas, S. (2021). The longitudinal trajectories of
online engagement over a full program. *Computers & Education, 175*,
104325. <https://doi.org/10.1016/j.compedu.2021.104325>

Tikka, S., López-Pernas, S., & Saqr, M. (2025). tna: An R package for
transition network analysis. *Applied Psychological Measurement*.
<https://doi.org/10.1177/01466216251348840>
