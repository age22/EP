#' Assembles a path and creates intermediary directories
#'
#' Given a character vector, it creates a non-OS specific path and makes all
#' intermediate directories necessary.
#'
#' @param parts character vector
#'
#' @return string
#' @export
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
