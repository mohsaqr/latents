# Plot a cross-level fit

`"profiles"` shows the individual profiles' means per indicator;
`"group_means"` the group classes' means of the group means;
`"composition"` the profile shares within each group class. Profiles and
classes are told apart by colour and point shape.

## Usage

``` r
# S3 method for class 'multilpa_cross_level'
plot(
  x,
  what = c("profiles", "group_means", "composition"),
  main = NULL,
  subtitle = NULL,
  ...
)
```

## Arguments

- x:

  A `multilpa_cross_level` fit.

- what:

  Which view.

- main, subtitle:

  Optional title and subtitle.

- ...:

  Unused.

## Value

A ggplot object. Raises `latents_missing_package` without ggplot2.
