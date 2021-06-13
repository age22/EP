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
#'
#' @examples
#' packages <- c("ggplot2", "xlsx")
#' load_packages(packages)
load_packages <-
  function(packages) {
    invisible(suppressPackageStartupMessages(
      lapply(packages, library, character.only = TRUE)))

    invisible(NULL)
  }

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
#'
#' @examples
#' possible_options <- c("EP", "verbose", "imputated", "ADNI")
#' assign_arguments(arguments = c("EP", "verbose"), options = possible_options)
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

#' Get the complementary allele sequence
#'
#' For a character vector of a genomic sequences it returns a character vector
#' with the complementary sequence.
#'
#' @param sequence character vector, with a nucleotide sequence.
#' @param mode string. If "dna" it complements "A" to "T". If "rna" it
#' complements "A" to "U".
#'
#' @return character vector
#' @export
#'
#' @examples
#' complement_alleles(c("A", "T", "G"))
complement_alleles <-
  function(sequence, mode = "dna") {

    # Creating complement function
    complement_function <- function(allele, .mode = mode) {
      if (allele == "G") {
        complement <- "C"
      }
      if (allele == "C") {
        complement <- "G"
      }
      if (allele == "A") {
        if (mode == "dna") {
          complement <- "T"
        } else if (mode == "rna") {
          complement <- "U"
        }
      }
      if (allele == "T" | allele == "U") {
        complement <- "A"
      }
      complement
    }
    # Using complement function on sequence elements
    sequence <- toupper(sequence)
    complementary <- sapply(sequence, complement_function)
    complementary <- unname(complementary)

    # Returning complementary sequence obtained
    complementary
  }

#' Assembles a path and creates intermediary directories
#'
#' Given a character vector, it creates a non-OS specific path and makes all
#' intermediate directories necessary.
#'
#' @param parts character vector
#'
#' @return string
#' @export
#'
#' @examples
#' create_paths(c("All", "the", "way", "to", "my", "directory"))
create_paths <-
  function(parts) {
    # Create path
    path <- do.call(file.path, as.list(parts))

    # If the intermediate folders doesn't exist, create them
    if (!dir.exists(file.path(path))) {
      dir.create(path, recursive = T)
    }
    path
  }

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
#'
#' @examples
#' rec_list(c(2,3))
rec_list <-
  function(x) {
    if (length(x) == 1) {
      vector("list", x)
    } else {
      lapply(1:x[1], function(...) rec_list(x[-1]))
    }
  }

#' Split string into a vector of components
#'
#' Wrapper around \code{\link[base]{strsplit}}, returning a vector instead of a list.
#'
#' @param string string
#' @param pattern string
#'
#' @return vector
#' @export
#'
#' @examples
#' strsplit2vector("hello_world", "_")
strsplit2vector <-
  function(string, pattern) {
    unlist(strsplit(x = string, split = pattern))
  }

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
      cat(system(command = command, intern = TRUE, ...), sep = "\n")
    }
    else if (is.windows) {
      cat(shell(cmd = command, intern = TRUE, translate = T, ...), sep = "\n")
    }
  }

