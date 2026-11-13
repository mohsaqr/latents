# Number one fit's profiles like another's

Profiles and latent classes have no natural order, so two fits of the
same model to different samples, waves or subgroups usually number them
differently. `align_labels()` renumbers the profiles of `fit`, and its
group classes in a two-level model, so that each takes the label of the
most similar profile (group class) of `reference`. Nothing is
re-estimated: the log likelihood, the estimates and every assignment are
the same, under new labels.

## Usage

``` r
align_labels(fit, reference)
```

## Arguments

- fit:

  A fit from
  [`lpa()`](https://pak.dynasite.org/latents/reference/lpa.md),
  [`lca()`](https://pak.dynasite.org/latents/reference/lca.md),
  [`multilpa()`](https://pak.dynasite.org/latents/reference/multilpa.md)
  or
  [`multilca()`](https://pak.dynasite.org/latents/reference/multilca.md),
  with or without membership covariates.

- reference:

  A fit of the same kind whose labels `fit` takes: the same indicators,
  coded the same way, and the same numbers of profiles and group
  classes.

## Value

`fit`, of the same class, with its profiles (and group classes)
renumbered: every table of
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
the assignments and the posterior probabilities follow the new labels.
Raises `latents_bad_argument` when either object is not a profile or
class model, or when the two fits differ in their indicators, how they
are coded (the message names each categorical or ordinal indicator whose
categories or their order differ), or their numbers of profiles or group
classes.

## Details

Profiles are matched on their measurement parameters (means, variances
and, when both fits estimate them, correlations; response probabilities;
ordinal and count parameters), measured in the reference's units, with
each profile's share as a tie-breaker. Every permutation is tried for up
to seven profiles, and each reference profile takes its nearest
unclaimed match above that. Group classes are then matched on the
profile mixture each carries. A membership-covariate fit may be aligned
to the same model without covariates: its share is the mean posterior
probability, and its membership logits are re-expressed against the new
last profile. This is the alignment the cluster bootstrap,
[`sensitivity()`](https://pak.dynasite.org/latents/reference/sensitivity.md)
and
[`pool_imputations()`](https://pak.dynasite.org/latents/reference/pool_imputations.md)
use.

## See also

[`compare_models()`](https://pak.dynasite.org/latents/reference/compare_models.md)
to compare fits on their criteria.

## Examples

``` r
first_days <- subset(student_esm, day <= 3)
last_days <- subset(student_esm, day > 3)
feelings <- c("happy", "exhausted")
early <- lpa(first_days, feelings, n_profiles = 2, n_starts = 2, seed = 1)
late <- lpa(last_days, feelings, n_profiles = 2, n_starts = 2, seed = 1)
# The later days' fit numbers the same two profiles the other way round.
get_results(early, "profiles")
#>   profile indicator     mean  variance standard_deviation mean_standard_error
#> 1       1     happy 5.565861 0.9550211          0.9772518          0.07763166
#> 2       1 exhausted 1.933072 0.6632697          0.8144137          0.06909314
#> 3       2     happy 4.077278 2.4600442          1.5684528          0.09053118
#> 4       2 exhausted 4.587211 2.0688916          1.4383642          0.11336925
#>   variance_standard_error
#> 1              0.10918578
#> 2              0.07583124
#> 3              0.18442118
#> 4              0.20080060
get_results(late, "profiles")
#>   profile indicator     mean  variance standard_deviation mean_standard_error
#> 1       1     happy 4.474893 2.5181753          1.5868760          0.04508037
#> 2       1 exhausted 4.022386 2.4217003          1.5561813          0.05167207
#> 3       2     happy 6.143550 0.6715874          0.8195043          0.04742078
#> 4       2 exhausted 1.322795 0.2250344          0.4743780          0.02741166
#>   variance_standard_error
#> 1              0.09858709
#> 2              0.10232709
#> 3              0.05341941
#> 4              0.01776771
get_results(align_labels(late, early), "profiles")
#>   profile indicator     mean  variance standard_deviation mean_standard_error
#> 1       1     happy 6.143550 0.6715874          0.8195043          0.04742078
#> 2       1 exhausted 1.322795 0.2250344          0.4743780          0.02741166
#> 3       2     happy 4.474893 2.5181753          1.5868760          0.04508037
#> 4       2 exhausted 4.022386 2.4217003          1.5561813          0.05167207
#>   variance_standard_error
#> 1              0.05341941
#> 2              0.01776771
#> 3              0.09858709
#> 4              0.10232709
```
