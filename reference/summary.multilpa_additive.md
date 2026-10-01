# Summarize an additive group-class fit

Summarize an additive group-class fit

## Usage

``` r
# S3 method for class 'multilpa_additive'
summary(object, level = 0.95, vcov_type = c("observed", "robust", "opg"), ...)

# S3 method for class 'summary_multilpa_additive'
print(x, digits = 4L, ...)

# S3 method for class 'summary_multilpa_additive'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- object:

  A `multilpa_additive` fit.

- level, vcov_type:

  As in
  [`get_results.multilpa_additive()`](https://pak.dynasite.org/latents/reference/get_results.multilpa_additive.md).

- ...:

  Unused.

- x:

  A `summary_multilpa_additive` object.

- digits:

  Significant digits.

- row.names, optional:

  Unused; part of the generic.

## Value

An object of class `summary_multilpa_additive`: a named list of the
`fit`, `parameters` and `group_classes` tables, printed by its method.
