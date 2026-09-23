# Number of independent units in a fitted latent transition model

Number of independent units in a fitted latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Ignored.

## Value

A single integer: the number of observed groups. Occasions within a
group are dependent by construction, so they are not independent
observations.

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
nobs(fit)
#> [1] 30
```
