# Plot information criteria across transition models

The criterion against the number of profiles, one line per number of
group classes (and structure), told apart by colour, shape and line
type.

## Usage

``` r
# S3 method for class 'latents_transition_enumeration'
plot(
  x,
  criterion = c("bic", "aic", "bic_individual"),
  main = NULL,
  subtitle = NULL,
  ...
)
```

## Arguments

- x:

  A `latents_transition_enumeration` result.

- criterion:

  `"bic"`, `"aic"` or `"bic_individual"`.

- main, subtitle:

  Optional title and subtitle.

- ...:

  Unused.

## Value

A ggplot object. Raises `latents_missing_package` without ggplot2.

## Examples

``` r
few_students <- subset(course_engagement, student <= 50)
grid <- enumerate_classes(few_students, c("browse", "lectures"),
                          "student", time = "sequence", n_profiles = 2:3,
                          n_group_classes = 1, n_starts = 2, seed = 1)
if (requireNamespace("ggplot2", quietly = TRUE)) plot(grid)
```
