# Compare fitted models on the same data

One row per model with its log likelihood, number of parameters,
information criteria, classification quality and convergence, so
alternative specifications (for example trajectory classes alone against
a growth mixture model) are compared without assembling tables by hand.
`delta_bic` is each model's BIC above the smallest, and `bic_weight` the
Schwarz weight, exp(-delta/2) normalized over the models compared
(Wagenmakers and Farrell 2004): the approximate posterior probability of
each model if one of them is true, with equal prior probabilities.

## Usage

``` r
compare_models(...)
```

## Arguments

- ...:

  Two or more fits from
  [`mixture_regression()`](https://pak.dynasite.org/latents/reference/mixture_regression.md)
  (with or without `random`). Name them to label the rows; unnamed fits
  are labelled by their expression.

## Value

A base `data.frame` of class `latents_comparison` with one row per
model: `model`, `n_classes`, `random`, `n_parameters`, `log_likelihood`,
`aic`, `bic`, `sabic`, `icl`, `entropy`, `smallest_share`, `delta_bic`,
`bic_weight`, `converged`. Raises `latents_bad_argument` for fits it
cannot compare and `latents_incomparable_models` for fits to different
data.

## Details

Information criteria are comparable only on the same data: every model
must be fitted to the same outcome, persons (or groups) and
observations, or the comparison is refused.

## References

Wagenmakers, E.-J., & Farrell, S. (2004). AIC model selection using
Akaike weights. *Psychonomic Bulletin & Review*, 11, 192–196.

## Examples

``` r
# Fixed trajectories against growth curves with a random intercept,
# for the first 60 students
few_students <- subset(growth_scores, student <= 60)
trajectories <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                                   id = "student", class_level = "group",
                                   seed = 1)
growth <- mixture_regression(score ~ wave, few_students, n_classes = 2,
                             id = "student", class_level = "group",
                             random = "intercept", seed = 1)
compare_models(trajectories = trajectories, growth = growth)
#> Model comparison
#> 
#> Model         Classes  Random     k    LogLik      BIC    ΔBIC  Weight  Entropy
#> ------------  -------  ---------  -  --------  -------  ------  ------  -------
#> trajectories        2  none       7  -1277.94  2584.53  196.11   0.000     0.94
#> growth              2  intercept  9  -1175.79  2388.43    0.00   1.000     0.97  <- best
```
