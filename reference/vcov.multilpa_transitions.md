# Standard errors are not available for a latent transition model

Standard errors are not available for a latent transition model

## Usage

``` r
# S3 method for class 'multilpa_transitions'
vcov(object, ...)

# S3 method for class 'multilpa_transitions'
parameter_inference(x, data = NULL, ...)

# S3 method for class 'multilpa_transitions'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A fitted `multilpa_transitions` model.
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) are base generics,
  so their first formal is `object`.

- ...:

  Ignored.

- x:

  The same fitted model, under the name this package's own verbs use;
  [`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
  is documented on this page too.

- data:

  Ignored; present for compatibility with the generic.

- parm:

  Ignored; present for compatibility with the generic.

- level:

  Ignored; present for compatibility with the generic.

## Value

Nothing; [`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`parameter_inference()`](https://pak.dynasite.org/latents/reference/parameter_inference.md)
all raise a `latents_no_inference` condition on a latent transition fit.
The analytic score and Jacobian this package uses for
observed-information and sandwich standard errors do not yet cover the
initial and transition multinomial logits, so no interval is reported
rather than an invalid one.
[`confint()`](https://rdrr.io/r/stats/confint.html) refuses explicitly
instead of letting
[`confint.default()`](https://rdrr.io/r/stats/confint.html) reach
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and refuse by accident.

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
tryCatch(confint(fit), latents_no_inference = function(condition) {
  conditionMessage(condition)
})
#> [1] "Standard errors are not available for a latent transition model. Read the estimates with get_results()."
```
