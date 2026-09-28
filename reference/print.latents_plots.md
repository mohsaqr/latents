# Print several plots

The value of `plot(x, what = "all")`, of
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
diagnostics result and of an enumeration plotted with `combine = FALSE`:
a list of ggplot objects, named by view. Printing draws each in turn.

## Usage

``` r
# S3 method for class 'latents_plots'
print(x, ...)
```

## Arguments

- x:

  A `latents_plots` list.

- ...:

  Ignored.

## Value

`x`, invisibly. Called for the side effect of drawing.
