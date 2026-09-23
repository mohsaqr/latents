# Print a class-enumeration grid

Print a class-enumeration grid

## Usage

``` r
# S3 method for class 'multilpa_enumeration'
print(x, ...)
```

## Arguments

- x:

  An `multilpa_enumeration` result.

- ...:

  Passed to the underlying `data.frame` printing.

## Value

The input, invisibly. Called for the side effect of printing the
candidate grid: one line per candidate with its class counts, covariance
model, log likelihood, parameter count, AIC, BIC under both sample-size
conventions, ICL counted over individuals, entropy at both levels, and
the diagnostics needed to trust a candidate (convergence, a bound
reached, and how many starts reached the best likelihood), then a count
of the candidates that did not converge.
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) give
every criterion.

## Examples

``` r
set.seed(7)
example_data <- data.frame(
  school = rep(seq_len(12), each = 10),
  score_a = rnorm(120), score_b = rnorm(120)
)
candidates <- enumerate_classes(
  example_data, c("score_a", "score_b"), "school", n_profiles = 1:2,
  n_group_classes = 1, n_starts = 2, seed = 1
)
print(candidates)
#> Class enumeration: 2 candidate models
#>  n_profiles n_group_classes model log_likelihood n_parameters      aic
#>           1               1   VVI      -327.6034            4 663.2068
#>           2               1   VVI      -323.0197            9 664.0394
#>  bic_groups bic_individual icl_individual profile_entropy group_entropy
#>    665.1464       674.3568       674.3568              NA            NA
#>    668.4036       689.1268       720.6548        0.810478            NA
#>  converged boundary n_best_replicated
#>       TRUE    FALSE                 2
#>       TRUE    FALSE                 2
#> No candidate is selected automatically. Compare one convention consistently.
```
