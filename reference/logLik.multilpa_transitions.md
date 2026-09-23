# Log likelihood of a fitted latent transition model

Log likelihood of a fitted latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
logLik(object, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.

- ...:

  Ignored.

## Value

A `logLik` object carrying the maximized observed-data log likelihood,
the free parameter count as `df`, and the number of groups as `nobs`.
Groups are the independent units, because a group's occasions are
dependent by construction in this model.

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
logLik(fit)
#> 'log Lik.' -418.2111 (df=11)
```
