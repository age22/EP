#' Gets the corresponding inheritance model
#'
#' This function gets a string with the corresponding inheritance model for each
#' SNP (dominant, recessive or additive) and returns the id of the SNP with the
#' first letter of the model (d, r or a) added to it. If the inheritance model
#' string is "unknown", it looks for the best model at the master_list object
#' and retrieves it.
#'
#' @param snps character vector, of snp id's
#' @param models character vector, with the inheritance model of the snps
#' @param master_list list, with all results from the analysis
#'
#' @returns a vector with all the snp_id + model
#'
#' @export
get_model <-
  function(snps, models, master_list) {
    vector_of_models <- vector(mode = "character", length = length(models))
    for (i in seq_along(snps)) {
      snp <- snps[i]
      model <- models[i]
      if (model == "unknown") {
        model <- names(master_list$All[[snp]]$Best_model)
      } else {
        if (model == "recessive") {
          model <- paste0(snp, "r")
        } else if (model == "dominant") {
          model <- paste0(snp, "d")
        } else if (model == "additive") {
          model <- paste0(snp, "a")
        }
      }
      vector_of_models[i] <- model
    }
    vector_of_models
  }


