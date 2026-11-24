# Plot every case as a line across the continuous indicators, one panel per profile

`plot_parallel()` plots each case as one line that joins its values on
the continuous indicators, with one panel per profile (a
parallel-coordinates plot).

## Usage

``` r
plot_parallel(x, main = NULL, subtitle = NULL)
```

## Arguments

- x:

  One of:

  - a profile fit from
    [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md) or
    [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md),
    with or without covariates, that has two or more continuous
    indicators and was fitted without centering or scaling;

  - a transition fit from
    [`lta()`](https://pak.dynasite.org/latents/reference/lta.md) with
    homogeneous transitions (class `multilpa_transitions`) that has two
    or more continuous indicators;

  - a
    [`diagnostics()`](https://pak.dynasite.org/latents/reference/diagnostics.md)
    result, which is drawn from the fit it checked.

- main, subtitle:

  Title and subtitle of the plot. `NULL` (default) keeps the view's own
  text.

## Value

A ggplot object.

Errors:

- `latents_nothing_to_plot` when the package returned `x` but this
  function has no plot of it: a fit with fewer than two continuous
  indicators, a centred or scaled fit, a general transition fit (class
  `multilpa_lta`), an additive or cross-level fit, a mixture regression,
  a growth mixture or an enumeration grid. For an object of the wrong
  kind, the message names the plot functions that apply to it;

- `latents_not_a_fit` when `x` is not a fit, grid or other result of the
  package, for example a data frame;

- `latents_missing_package` when ggplot2 is not installed.

## Details

The x-axis lists the continuous indicators and the y-axis is the
standardized value (z), so that indicators on different scales can share
one axis. Each panel holds the cases whose most likely profile it is,
drawn in the profile's colour; a line is fainter the less certain the
case's assignment. A case with a missing value has its line broken at
the gap, and a value with no answered neighbour is drawn as a point. The
cases of the other profiles are drawn in grey behind. The dark line with
markers is the mean of the panel's profile.

The view needs two or more continuous indicators and draws the stored
data of the fit. A fit with person or grand centering or scaling is
refused, because its stored data are not on the scale of the model. For
a transition fit, a case is one answered occasion of a unit, counted for
its most likely state.

## See also

[`plot_raincloud()`](https://pak.dynasite.org/latents/reference/plot_raincloud.md)
and
[`plot_pairs()`](https://pak.dynasite.org/latents/reference/plot_pairs.md)
for the other views of the cases;
[`plot_profiles()`](https://pak.dynasite.org/latents/reference/plot_profiles.md)
for the profile means.

## Examples

``` r
fit <- lpa(faithful, c("eruptions", "waiting"), n_profiles = 2,
           n_starts = 2, seed = 1)
plot_parallel(fit)
```
