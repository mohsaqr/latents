# Conditions raised by this package never carry data in their call.
#
# A fitting verb reached through do.call() records its arguments by value, so
# an error raised while validating the input carried the whole data frame in
# `conditionCall()`; R's own printing truncates the call, but a host that
# deparses it in full streams every value (review 2026-10-27). Entry points
# therefore evaluate calls whose arguments are symbols, and every condition
# that is stored on a result or raised again is bounded first.

#' Bind the arguments of a call to symbols in an environment
#'
#' Unless `values` are supplied, each argument is evaluated once in `evaluate_in` (the caller's frame, so
#' an expression's side effects land where a promise's would) and its value
#' is bound in `env`: under the caller's own name when the argument is a
#' symbol, otherwise under `.<name>`, with a unique suffix for collisions.
#' The returned call therefore records
#' symbols only, never the data, while evaluating to the same values.
#' @param call A call.
#' @param env The environment the call will be evaluated in.
#' @param evaluate_in The environment the caller's expressions belong to.
#' @param values Optional already-evaluated arguments, in call order. Wrappers
#'   that have forced their promises supply these instead of evaluating the
#'   caller's expressions again.
#' @return The call with every argument replaced by a symbol bound in `env`.
#' @noRd
.latents_bind_values <- function(call, env, evaluate_in, values = NULL) {
  stopifnot("`call` must be a call" = is.call(call),
            "`env` must be an environment" = is.environment(env),
            "`evaluate_in` must be an environment" = is.environment(evaluate_in),
            "`values` must contain one value per argument" =
              is.null(values) || (is.list(values) && length(values) == length(call) - 1L))
  indices <- seq_along(call)[-1L]
  names <- names(call) %||% rep("", length(call))
  # Generated names must not shadow a later caller symbol or the function
  # already bound in env. Repeated symbols get separate slots too.
  symbols <- vapply(indices, function(index) {
    if (is.symbol(call[[index]])) as.character(call[[index]]) else ""
  }, character(1))
  reserved <- ls(env, all.names = TRUE)
  invisible(lapply(indices, function(index) {
    argument <- call[[index]]
    value <- if (!is.null(values)) values[[index - 1L]] else
      if (is.language(argument)) eval(argument, evaluate_in) else argument
    name <- if (is.symbol(argument)) as.character(argument) else
      if (nzchar(names[[index]])) paste0(".", names[[index]]) else
        sprintf(".argument_%d", index)
    occupied <- unique(c(reserved, if (!is.symbol(argument)) symbols))
    if (name %in% occupied) {
      name <- utils::tail(make.unique(c(occupied, name)), 1L)
      # A suffix can itself be a caller's symbol reserved for a later arg.
      name <- utils::tail(make.unique(c(unique(c(reserved, symbols)), name)), 1L)
    }
    reserved <<- c(reserved, name)
    assign(name, value, envir = env)
    call[[index]] <<- as.symbol(name)
  }))
  call
}

#' A call whose deparse is bounded, for storing or raising again
#'
#' A function object at the head becomes `<function>`; a value argument whose
#' one-line deparse is longer than `max_chars` becomes a placeholder naming
#' its class and size, such as `<data.frame[400000 x 9]>`. Symbols and
#' expressions are kept as the caller wrote them.
#' @param call A call, or anything else (returned unchanged).
#' @param max_chars The longest value kept verbatim.
#' @return The bounded call.
#' @noRd
.latents_bounded_call <- function(call, max_chars = 80L) {
  if (!is.call(call)) return(call)
  parts <- as.list(call)
  if (is.function(parts[[1L]])) parts[[1L]] <- as.symbol("<function>")
  rest <- seq_along(parts)[-1L]
  parts[rest] <- lapply(parts[rest], function(value) {
    if (is.language(value) || is.null(value)) return(value)
    text <- deparse(value, width.cutoff = 500L, nlines = 2L)
    if (length(text) == 1L && nchar(text) <= max_chars) return(value)
    as.symbol(.latents_value_placeholder(value))
  })
  as.call(parts)
}

.latents_value_placeholder <- function(value) {
  dims <- dim(value)
  size <- if (is.null(dims)) format(length(value)) else paste(dims, collapse = " x ")
  sprintf("<%s[%s]>", class(value)[1L], size)
}

#' A condition with a bounded call
#' @param condition A condition.
#' @return The same condition with `.latents_bounded_call()` applied.
#' @noRd
.latents_bounded_condition <- function(condition) {
  stopifnot("`condition` must be a condition" = inherits(condition, "condition"))
  if (is.call(conditionCall(condition))) {
    condition$call <- .latents_bounded_call(conditionCall(condition))
  }
  condition
}
