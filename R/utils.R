# Internal helpers. No Rd page: nothing here is exported.

# Polyfill `%||%` for R < 4.4. Base R 4.4 added it; this package declares
# R (>= 4.1.0), so the base version cannot be relied on. Defined
# unconditionally because the package namespace is searched first, and base's
# version on R >= 4.4 is functionally identical, so this is a harmless shadow.
`%||%` <- function(a, b) if (is.null(a)) b else a
