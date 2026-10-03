# Print a latents result table

Prints the columns people report, with readable labels and numbers; the
object itself is unchanged, a base data frame with every numeric column.

## Usage

``` r
# S3 method for class 'latents_table'
print(x, n = 20L, ...)

# S3 method for class 'latents_table'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `latents_table`.

- n:

  The most rows printed.

- ...:

  Unused.

- row.names, optional:

  Passed to the data frame method.

## Value

`x`, invisibly.

For [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), the
same table as a plain base data frame.
