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


