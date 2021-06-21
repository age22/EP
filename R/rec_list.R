#' Create empty recursive list from a numerical vector
#'
#' This function creates and empty list for a numerical vector so that every
#' position corresponds to the number of elements inside that specific recursion
#' level being \code{length(vector)} the number of recursions levels of the list.
#'
#' @param x numerical vector
#'
#' @return list
#' @export
rec_list <-
  function(x) {
    if (length(x) == 1) {
      vector("list", x)
    } else {
      lapply(1:x[1], function(...) rec_list(x[-1]))
    }
  }
