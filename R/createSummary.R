#' Create Summary
#'
#' @param variable a variable
#' @param group a group
#'
#' @return a summary table
#' @export
createSummary <- function(variable, group) {
  capture.output(psych::describeBy(variable, group, digits = 2))
}

