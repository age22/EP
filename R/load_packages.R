#' Load and attach libraries from vector of names.
#'
#' Convenience function that loads and attaches all libraries named in a
#' character vector silently.
#'
#' @param packages character vector, with the names of the libraries to be
#' loaded
#'
#' @return invisible NULL
#' @export
load_packages <-
  function(packages) {
    invisible(suppressPackageStartupMessages(
      lapply(packages, library, character.only = TRUE)))

    invisible(NULL)
  }
