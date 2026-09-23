# Print a fitted latent transition model

Reports the model's size, how the occasions are laid out, its fit and
its convergence, and names the verbs that return the fitted quantities.

## Usage

``` r
# S3 method for class 'multilpa_transitions'
print(x, rows = 20L, ...)
```

## Arguments

- x:

  A fitted `multilpa_transitions` model.

- rows:

  How many rows of the printed table to show before truncating.

- ...:

  Ignored.

## Value

`x`, invisibly. Called for the side effect of printing the class counts,
the occasion layout, the log likelihood with the information criteria,
the convergence and restart diagnostics, and the verbs that return the
fitted quantities.

## See also

[`get_results()`](https://pak.dynasite.org/latents/reference/get_results.md),
[`lta()`](https://pak.dynasite.org/latents/reference/lta.md).

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
print(fit)
#> Latent transition model: 2 profiles, 1 group class
#> 150 observations in 30 groups, up to 5 occasions (balanced)
#> Log likelihood: -418.211098; 11 parameters; BIC (groups): 873.8354
#> Converged: TRUE after 116 iterations; best of 2 starts
#> 
#>  profile    score_a     score_b     count proportion
#>        1 0.09137155 -0.09015126 135.99729 0.90664862
#>        2 0.76447843  0.89468659  14.00271 0.09335138
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
