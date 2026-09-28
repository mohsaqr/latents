# Print a fitted multilevel latent profile model

Print a fitted multilevel latent profile model

## Usage

``` r
# S3 method for class 'multilpa'
print(x, rows = 20L, ...)
```

## Arguments

- x:

  An `multilpa` model.

- rows:

  How many rows of the printed table to show before truncating.

- ...:

  Reserved for compatibility with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The input model, invisibly. Called for the side effect of printing the
class counts, the sample sizes and covariance specification, the log
likelihood with AIC and group-level BIC, the convergence and restart
diagnostics, and any blocks the fit held fixed.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
fit <- multilpa(example_data, c("score_a", "score_b"), "school",
                n_profiles = 2, n_group_classes = 1, n_starts = 2, seed = 1)
print(fit)
#> Two-level latent profile analysis: 2 profiles, 1 group class
#> 120 individuals in 12 groups; varying diagonal residual covariance (VVI)
#> Log likelihood: -323.019688 | AIC: 664.039 | BIC (groups): 668.404
#> Converged: TRUE | iterations: 30 | best start: 2/2
#> 
#>  profile    score_a     score_b     count proportion
#>        1  0.3163698 -0.05291982 106.38528   0.886544
#>        2 -1.1135958  0.79091589  13.61472   0.113456
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
