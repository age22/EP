#' Run programs in the shell
#'
#' Wrapper around the system() and shell() functions that works on Rmarkdown
#'
#' @param command character string, quoted command
#' @param ... other arguments accepted by the system() function apart from
#' intern, which is always set to TRUE.
#'
#' @return invisible NULL
#' @export
run_shell <-
  function(command, ...) {
    is.windows <- ifelse(.Platform$OS.type == "windows", TRUE, FALSE)
    is.unix <- ifelse(.Platform$OS.type == "unix", TRUE, FALSE)
    if (is.unix) {
      cat(system(command = command, ...), sep = "\n")
    }
    else if (is.windows) {
      cat(shell(cmd = command, translate = T, ...), sep = "\n")
    }
  }

