# Fitted parameters of a latent transition model

Fitted parameters of a latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
coef(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Reserved for compatibility with
  [`coef()`](https://rdrr.io/r/stats/coef.html).

## Value

A named numeric vector of every free parameter on its natural scale:
profile means and variances, categorical response probabilities where
present, each group class's initial profile probabilities, its
transition probabilities, and the group-class probabilities.

Names follow the package-wide `level.parameter.outcome.term` grammar
shared with
[`coef.multilpa()`](https://pak.dynasite.org/latents/reference/coef.multilpa.md),
so one pattern matches across fit classes: for example
`measurement.mean.profile_1.reading` and
`profile.transition_probability.profile_2.group_class_1:profile_1`,
whose term names the origin the move is from. `level`, `parameter` and
`outcome` never contain a dot, so everything after the third dot is the
term and the name parses back into the columns
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
reports even when an indicator name itself contains a dot.

No standard errors accompany these:
[`vcov.multilpa_transitions()`](https://pak.dynasite.org/latents/reference/vcov.multilpa_transitions.md)
refuses rather than returning an invalid matrix.
[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md)
gives the same quantities as tidy tables, which is the form to prefer.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  person = rep(seq_len(30), each = 5), wave = rep(seq_len(5), times = 30)
)
example_data$score_a <- stats::rnorm(nrow(example_data))
example_data$score_b <- stats::rnorm(nrow(example_data))
fit <- lta(example_data, c("score_a", "score_b"), "person",
                       n_profiles = 2, time = "wave", n_starts = 2, seed = 1)
coef(fit)
#>                               measurement.mean.profile_1.score_a 
#>                                                     9.137155e-02 
#>                               measurement.mean.profile_1.score_b 
#>                                                    -9.015126e-02 
#>                               measurement.mean.profile_2.score_a 
#>                                                     7.644784e-01 
#>                               measurement.mean.profile_2.score_b 
#>                                                     8.946866e-01 
#>                           measurement.variance.profile_1.score_a 
#>                                                     7.803627e-01 
#>                           measurement.variance.profile_1.score_b 
#>                                                     1.086468e+00 
#>                           measurement.variance.profile_2.score_a 
#>                                                     1.159944e+00 
#>                           measurement.variance.profile_2.score_b 
#>                                                     5.082751e-01 
#>              profile.initial_probability.profile_1.group_class_1 
#>                                                     9.066520e-01 
#>              profile.initial_probability.profile_2.group_class_1 
#>                                                     9.334804e-02 
#> profile.transition_probability.profile_1.group_class_1:profile_1 
#>                                                     9.999980e-01 
#> profile.transition_probability.profile_2.group_class_1:profile_1 
#>                                                     2.012626e-06 
#> profile.transition_probability.profile_1.group_class_1:profile_2 
#>                                                     1.000000e-10 
#> profile.transition_probability.profile_2.group_class_1:profile_2 
#>                                                     1.000000e+00 
#>                                  group.probability.group_class_1 
#>                                                     1.000000e+00 
```
