#' Store interactions of SNP-SNP
#'
#' @param snp_dataset snp dataframe
#' @param master_list master list with all results
#'
#' @return data.frame, the master_list with all updated results
#' @export
store_interactions <-
  function(snp_dataset, master_list) {

    # Getting all models of a certain SNP1
    snps1 <- unique(snp_dataset$snp1)
    models1 <- unique(snp_dataset$model1)
    models_per_snp1 <- lapply(snps1, function(n) grep(pattern = n, models1, value = T))
    names(models_per_snp1) <- snps1
    interactions_list <- models_per_snp1

    # Adding all SNP2 that interact with a certain SNP1 model
    snp2_per_model <- lapply(unlist(interactions_list), function(x) unique(snp_dataset$snp2[snp_dataset$model1 == x]))
    names(snp2_per_model) <- unlist(interactions_list)

    for (snp1 in names(interactions_list)) {
      models1_list <- vector("list", length = length(interactions_list[[snp1]]))
      models1 <- interactions_list[[snp1]]
      for (i in seq_along(models1)) {
        model1 <- models1[i]
        names(models1_list)[i] <- model1
        models1_list[[i]] <- snp2_per_model[[model1]]
      }
      interactions_list[[snp1]] <- models1_list
    }

    ## Adding all models of a certain SNP2 that interact with a given SNP1 model
    for (snp1 in names(interactions_list)) {
      for (model1 in names(interactions_list[[snp1]])) {
        snps2_list <- vector("list", length = length(interactions_list[[snp1]][[model1]]))
        snps_2 <- interactions_list[[snp1]][[model1]]
        for (i in seq_along(snps_2)) {
          snp2 <- snps_2[i]
          names(snps2_list)[i] <- snp2
          snp2_models <- snp_dataset$model2[snp_dataset$model1 == model1 & snp_dataset$snp2 == snp2]
          snps2_list[[i]] <- snp2_models
        }
        interactions_list[[snp1]][[model1]] <- snps2_list
      }
    }

    for (snp1 in names(interactions_list)) {
      for (model1 in names(interactions_list[[snp1]])) {
        for (snp2 in names(interactions_list[[snp1]][[model1]])) {
          models2_list <- vector("list", length = length(interactions_list[[snp1]][[model1]][[snp2]]))
          models_2 <- interactions_list[[snp1]][[model1]][[snp2]]
          for (i in seq_along(models_2)) {
            model2 <- models_2[i]
            names(models2_list)[i] <- model2
            models2_list[[i]] <- list("Name" = character() , "Summary" = character(), "SF" = numeric(), "Significant" = logical())
          }
          interactions_list[[snp1]][[model1]][[snp2]] <- models2_list
        }
      }
    }

    # Creating the containers that will store the interactions from the analysis

    for (snp in names(interactions_list)) {
      for (dataset in names(master_list)) {
        master_list[[dataset]][[snp]]$Interactions$SNP$Reported <- interactions_list[[snp]]
      }
    }
    master_list
  }
