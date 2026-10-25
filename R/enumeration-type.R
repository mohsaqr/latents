#' Validate an explicit enumeration type and its routing arguments
#' @param type NULL or a single model type.
#' @param time Occasion column, or NULL.
#' @param family Enumeration family.
#' @return NULL or the matched model type.
#' @noRd
.latents_check_enumeration_type <- function(type, time, family) {
  if (is.null(type)) return(NULL)
  invalid <- function() {
    stop(errorCondition(paste(
      "`type` must be NULL or a single string: \"lca\", \"lpa\" or \"lta\".",
      "Choose one type, or leave `type = NULL` for the existing routing."),
      class = "latents_bad_argument", call = NULL))
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) invalid()
  type <- tryCatch(match.arg(type, c("lca", "lpa", "lta")),
                    error = function(error) invalid())
  if (!identical(family, "profiles")) {
    stop(errorCondition(sprintf(paste(
      "`type = \"%s\"` conflicts with a non-\"profiles\" `family`.",
      "Set `family = \"profiles\"`, or leave `type = NULL` for a group-class family grid."),
      type), class = "latents_bad_argument", call = NULL))
  }
  if (type == "lta" && is.null(time)) {
    stop(errorCondition(paste(
      "`type = \"lta\"` requires `time`.",
      "Supply the occasion column with `time`, or choose `type = \"lca\"` or \"lpa\"."),
      class = "latents_bad_argument", call = NULL))
  }
  if (type != "lta" && !is.null(time)) {
    stop(errorCondition(sprintf(paste(
      "`type = \"%s\"` conflicts with `time`.",
      "Leave `time = NULL`, or use `type = \"lta\"` for transitions."), type),
      class = "latents_bad_argument", call = NULL))
  }
  type
}

#' Infer a profile grid's measurement label without changing its arguments
#' @param vars Indicator names.
#' @param extra Candidate arguments.
#' @return lca, lpa or mixed.
#' @noRd
.latents_enumeration_measurement_type <- function(vars, extra) {
  non_gaussian <- unlist(lapply(extra[c("categorical", "ordinal", "count")],
    function(indicators) if (is.character(indicators)) indicators else character()), use.names = FALSE)
  if (length(vars) > 0L && all(vars %in% non_gaussian)) "lca" else
    if (any(vars %in% non_gaussian)) "mixed" else "lpa"
}

#' Show the resolved model type in an enumeration header
#' @param x An enumeration, including objects saved before type was recorded.
#' @return A header suffix, or an empty string for older objects.
#' @noRd
.latents_enumeration_type_label <- function(x) {
  if (is.null(x$type)) "" else sprintf(" (%s)", toupper(x$type))
}
