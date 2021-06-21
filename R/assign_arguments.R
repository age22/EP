#' Assigns input arguments
#'
#' From a list of available options, assigns all given arguments as TRUE
#' and the rest as FALSE, and gives error if any of the inputted arguments
#' doesn't match the possible ones.
#'
#' @param arguments character vector, with inputted arguments names.
#' @param options character vector, with all possible options.
#'
#' @return invisible NULL
#' @export
assign_arguments <-
  function(arguments, options) {
    # Takes arguments inputted, compares it with possible options and sets
    # possible options inputted to TRUE and the remaning options as FALSE
    remaining_options <- options
    for (argument in arguments) {
      # Compare each inputted argument with possible options
      if (!argument %in% options) {
        error <- c("argument %s not one of the available script options.",
                   "Possible options are %s, %s and %s.")
        error <- paste(error, collapse = " ")
        stop(sprintf(error, argument, options[1], options[2], options[3]))
      }

      remaining_options <- setdiff(remaining_options, argument)
      argument <- sub("--", "", argument) # Remove double dash from name
      assign(argument, value = TRUE, envir = .GlobalEnv)
    }
    # Set remaining options as FALSE (default behaviour)
    for (option in remaining_options) {
      option <- sub("--", "", option)
      assign(option, FALSE, envir = .GlobalEnv)
    }
    invisible(NULL)
  }
