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
fit <- lta(subset(course_engagement, student <= 40),
           c("browse", "lectures", "forum_read"), "student",
           n_profiles = 2, time = "sequence", n_starts = 2, seed = 1)
print(fit)
#> Latent transition model: 2 profiles, 1 group class
#> 547 observations in 40 groups, up to 15 occasions (unbalanced, observed grid)
#> Log likelihood: -2025.406100; 15 parameters; BIC (groups): 4106.1454
#> Converged: TRUE after 8 iterations; best of 2 starts
#> 
#>  profile     browse   lectures forum_read    count proportion
#>        1  0.5552352  0.4372193  0.6089409 359.9133  0.6579768
#>        2 -0.8706168 -0.6264837 -0.9319730 187.0867  0.3420232
#> 
#> Variances and standard errors: get_results(x, "profiles"). 
#> Every other table: get_results(x, what = ), or get_results(x, "all").
```
