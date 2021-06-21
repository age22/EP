#' Split string into a vector of components
#'
#' Wrapper around \code{\link[base]{strsplit}}, returning a vector instead of a list.
#'
#' @param string string
#' @param pattern string
#'
#' @return vector
#' @export
strsplit2vector <-
  function(string, pattern) {
    unlist(strsplit(x = string, split = pattern))
  }
