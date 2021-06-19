
# PREPARATION -------------------------------------------------------------

#' Generate inheritances for a given set of snps in a dataset.
#'
#' Given a SNP_set object and a dataframe with the genotypes for each case for
#' those snps, recodify in new variables the genotypes into each one of three
#' inheritance models (additive, dominant and recessive.)
#'
#' @param snps SNP_set object
#' @param data dataframe, where \code{\link[base]{names}(snps)} corresponds to
#' variables in data.
#'
#' @return The same dataframe as in data but with 3 * \code{length(snps)} new columns added to it. Once for each inheritance model per snp.
#' @export
generate_inheritances <-
  function(snps, data) {
    for (snp in names(snps)) {
      # snp_a, snp_d or snp_r depending on the mode of inheritance.
      snp_recessive <- paste0(snp, "r")
      snp_dominant <- paste0(snp, "d")
      snp_additive <- paste0(snp, "a")

      # Getting all dataset genotypes for a particular snp
      genotype <- data[[snp]]

      # Retrieving snp object from the list of objects and its major and minor allele.
      SNP_object <- snps[[snp]]
      major_allele <- SNP_object$major_allele
      minor_allele <- SNP_object$minor_allele

      # Conditions to check if a genotype has a determined allele
      has_major_allele <- grepl(major_allele, genotype)
      has_minor_allele <- grepl(minor_allele, genotype)

      # Getting homozygotes genotypes
      homoz_major <- paste0(major_allele, major_allele)
      homoz_minor <- paste0(minor_allele, minor_allele)


      data[[snp_recessive]] <- ifelse(is.na(genotype), NA, ifelse(has_major_allele, 0, 1))
      data[[snp_dominant]] <- ifelse(is.na(genotype), NA, ifelse(has_minor_allele, 1, 0))
      data[[snp_additive]] <- ifelse(is.na(genotype), NA, ifelse(genotype == homoz_major, 0, ifelse(genotype == homoz_minor, 2, 1)))
    }
    data
  }



# ANALYSIS ----------------------------------------------------------------


perform_subset_analysis <-
  function(snp, dataset, covariate, covariates, verbose) {

    #Preparing input for the glm
    snp_model <- names(master_list[["All"]][[snp]][["Best_model"]])
    gene <- master_list[["All"]][[snp]][["Gene"]]
    covariate <- sub("_.*$", "", dataset)
    SF <- master_list[["All"]][[snp]][["Interactions"]][["Other_covariates"]][[covariate]][["SF"]]
    predictors <- c(covariates, snp_model)
    formula <- as.formula(paste("Diag", paste(predictors, collapse = " + "), sep = " ~ "))

    # Defining arguments to pass to the glm function
    args <-  list(formula = formula, family = quasibinomial("logit"), data = as.name(dataset))

    # Execute glm function
    linear_model <- do.call(glm, args)
    OR_summary <- Coeff.OR2(linear_model)[, c(3,6)]
    if (verbose) {
      cat("---------------------------------\n")
      print(gene); print(snp_model); print(dataset); print(SF)
      if (grepl("a", snp_model, fixed = T)) {
        print(tail(OR_summary, n = 2))
      } else {
        print(tail(OR_summary, n = 1))
      }
    }
  }

#' @export
print_ME <-
  function(list, corrected = F) {

    if (!corrected) {
    # Getting just the significant Main Effects (ME) p-value column
    sig_ME <- lapply(list$All, function(snp_list) {
      best_model <- names(snp_list[["Best_model"]]);
      main_effects <- snp_list[["Main_effects"]][[best_model]];
      main_effects[main_effects[,6] < 0.05,][,6]})

    # Remove snp info and just keep the pvalues for each significant term
    terms <- unname(sig_ME)

    } else {

      # Getting all Main Effects (ME) p-value column
      sig_ME <- lapply(list$All, function(snp_list) {
        lapply(snp_list$Main_effects, function(main_effects_per_model) {
          main_effects_per_model[,6]
        })})

      # Remove snp info and just keep the pvalues for each term
      terms <- lapply(unname(sig_ME), function(x) {unname(x)})
    }

    terms <- unlist(terms, use.names = T)

    # Keep the pvalues and the term separately
    p_values <- unname(terms)
    terms <- names(terms)

    # Keep only the snp_model terms
    snps <- grep("rs[0-9]*[a-z][0-2]$", terms, value = T)
    snp_indexes <- grep("rs[0-9]*[a-z][0-2]$", terms, value = F)

    # Merge again the pvalues and its corresponding term.
    p_values <- p_values[snp_indexes]

    if (corrected) {
      p_values <- p.adjust(p_values, method = "bonferroni")
      significant <- which(p_values < 0.05)
      p_values <- p_values[significant]
      snps <- snps[significant]
      is_best_model <- sapply(snps, function(snp_model_lvl) {
                              snp_model <- sub("[0-2]$", "", snp_model_lvl)
                              snp <- sub("[a-z]$", "", snp_model)
                              best_model <- list$All[[snp]]$Best_model
                              best_model <- names(best_model)
                              identical(snp_model, best_model)
      })

      p_values <- p_values[is_best_model]
      snps <- snps[is_best_model]
    }

    names(p_values) <- snps

    cat("\n")
    print(round(p_values, 3))
    cat("\n")

    # For each of these snps find the gene they belong to and count them
    sig_genes <- vector(mode = "character", length = length(snps))
    for (n in seq_along(snps)) {
      snp_model <- sub("[0-2]$", "", snps[n])
      snp <- sub("[a-z]$", "", snp_model)
      gene <- list_of_objects[[snp]][["gene"]]
      sig_genes[n] <- gene
    }
    cat("\n")
    print(sort(table(sig_genes), decreasing = T))
    cat("\n")

    invisible(NULL)
}

#' @export
print_INT <-
  function(list, covariates, corrected = F) {

    # Getting just the significant Interaction (INT) p-value column
    if (!corrected) {

    sig_INT <- lapply(list$All, function(snp_list) {
      interactions <- snp_list$Interactions$Other_covariates
      lapply(interactions, function(int_list){
        interaction <- int_list$Summary;
        interaction[interaction[,6] < 0.05,][,6]})})

    # Remove snp info and just keep the pvalues for each significant term
    terms <- lapply(unname(sig_INT), function(x) {unname(x)})


    # Getting all Interaction (INT) p-value column
    } else {

      sig_INT <- lapply(list$All, function(snp_list) {
        interactions <- snp_list$Interactions$Other_covariates
        lapply(interactions, function(int_list){
          interaction <- int_list$Summary;
          interaction[,6]})})

    # Remove snp info and just keep the pvalues for each term
      terms <- lapply(unname(sig_INT), function(x) {unname(x)})
    }

    terms <- unlist(terms, use.names = T)

    # Keep the pvalues and the term separately
    p_values <- unname(terms)
    terms <- names(terms)

    # Keep only the interaction terms
    snps <- grep("rs[0-9]*[a-z][0-2]:", terms, value = T)
    snp_indexes <- grep("rs[0-9]*[a-z][0-2]:", terms, value = F)

    # Merge again the pvalues and its corresponding term.
    p_values <- p_values[snp_indexes]

    if (corrected) {
      p_values <- p.adjust(p_values, method = "bonferroni")
      significant <- which(p_values < 0.05)
      p_values <- p_values[significant]
      snps <- snps[significant]
    }

    names(p_values) <- snps
    if (length(p_values > 0)) {
      cat("\n")
      print(round(p_values, 3))
      cat("\n")

      # For each of the covariate-snps interaction find the gene they belong to and count them
      for (covariate in covariates) {
        snp_cov <- grep(covariate, snps, value = T)
        sig_genes <- vector(mode = "character", length = length(snp_cov))
        for (n in seq_along(snp_cov)) {
          snp_model <- sub(":.*$", "", snp_cov[n])
          snp <- sub("[a-z][0-2]$", "", snp_model)
          gene <- list_of_objects[[snp]][["gene"]]
          sig_genes[n] <- gene
        }
        cat("\n")
        print(covariate)
        print(sort(table(sig_genes), decreasing = T))
        cat("\n")
      }
    } else {
      cat("No significant interactions found")
    }

    invisible(NULL)
  }

#' @export
print_SNP_SNP <-
  function(master_list, list_of_snp_pairs, list_of_APOE, corrected = F) {

    # Getting just the significant SNP-SNP interaction (SNP) p-value column
    if (!corrected) {
      sig_SNP <- lapply(list_of_snp_pairs, function(snp_pair) {
        snp1_model <- snp_pair[1]
        snp1 <- sub(snp1_model, pattern = "[a-z]$", replacement = "")
        snp2_model <- snp_pair[2]
        snp2 <- sub(snp2_model, pattern = "[a-z]$", replacement = "")
        interaction <- master_list$All[[snp1]][["Interactions"]][["SNP"]][["Reported"]][[snp1_model]][[snp2]][[snp2_model]][["Summary"]]
        if (!is.null(interaction[[1]])) {
          result <- interaction[,6][interaction[,6] < 0.05]
        } else {
          result <- NULL
        }
        result
      })
      sig_APOE <- lapply(list_of_APOE, function(snp_pair) {
        snp1_model <- snp_pair[1]
        snp1 <- sub(snp1_model, pattern = "[a-z]$", replacement = "")
        covariate <- snp_pair[2]
        interaction <- master_list$All[[snp1]][["Interactions"]][["Other_covariates"]][[covariate]][["Summary"]]
        if (!is.null(interaction[[1]])) {
          result <- interaction[,6][interaction[,6] < 0.05]
        } else {
          result <- NULL
        }
        result
      })


    # Getting all SNP-SNP interaction (SNP) p-value column
    } else {
      sig_SNP <- lapply(list_of_snp_pairs, function(snp_pair) {
        snp1_model <- snp_pair[1]
        snp1 <- sub(snp1_model, pattern = "[a-z]$", replacement = "")
        snp2_model <- snp_pair[2]
        snp2 <- sub(snp2_model, pattern = "[a-z]$", replacement = "")
        interaction <- master_list$All[[snp1]][["Interactions"]][["SNP"]][["Reported"]][[snp1_model]][[snp2]][[snp2_model]][["Summary"]]
        if (!is.null(interaction[[1]])) {
          result <- interaction[,6]
        } else {
          result <- NULL
        }
        result
      })
      sig_APOE <- lapply(list_of_APOE, function(snp_pair) {
        snp1_model <- snp_pair[1]
        snp1 <- sub(snp1_model, pattern = "[a-z]$", replacement = "")
        covariate <- snp_pair[2]
        interaction <-  master_list$All[[snp1]][["Interactions"]][["Other_covariates"]][[covariate]][["Summary"]]
        if (!is.null(interaction[[1]])) {
          result <- interaction[,6]
        } else {
          result <- NULL
        }
        result
      })
    }

    terms <- unlist(sig_SNP, use.names = T)
    terms2 <- unlist(sig_APOE, use.names = T)
    terms <- c(terms, terms2)

    # Keep the pvalues and the term separately
    p_values <- unname(terms)
    terms <- names(terms)

    # Keep only the interaction terms
    snps <- grep("rs[0-9]*[a-z][0-2]:", terms, value = T)
    snp_indexes <- grep("rs[0-9]*[a-z][0-2]:", terms, value = F)

    # Merge again the pvalues and its corresponding term.
    p_values <- p_values[snp_indexes]

    if (corrected) {
      p_values <- p.adjust(p_values, method = "bonferroni")
      significant <- which(p_values < 0.05)
      p_values <- p_values[significant]
      snps <- snps[significant]
    }

    names(p_values) <- snps

    if (length(p_values > 0)) {
      cat("\n")
      print(round(p_values, 3))
      cat("\n")

      # For each of the SNP-SNP interaction find the gene they belong to and count them
      sig_genes <- vector(mode = "character", length = length(snps))
      for (n in seq_along(snps)) {
        snp_model1 <- sub(":.*$", "", snps[n])
        snp1 <- sub("[a-z][0-2]$", "", snp_model1)
        snp_model2 <- sub("^.*:", "", snps[n])
        snp2 <- sub("[a-z][0-2]$", "", snp_model2)
        gene1 <- list_of_objects[[snp1]][["gene"]]
        gene2 <- list_of_objects[[snp2]][["gene"]]
        if (snp_model2 == "E4statusE4+") {
          gene2 <- "APOE"
        }
        sig_genes[n] <- paste(gene1, gene2, sep = "-")
      }
      cat("\n")
      print(sort(table(sig_genes), decreasing = T))
      cat("\n")

    } else {
      cat("No significant interactions found")
    }
  }



store_possible_interactions <-
  function(.master_list, pairs) {
    for (n in seq_along(.master_list)) {
      for (i in seq_along(table(sapply(pairs, function(x) x[1])))) {
        snp_length <- unname(table(sapply(pairs, function(x) x[1]))[i])
        snp_name <- names(table(sapply(pairs, function(x) x[1]))[i])
        .master_list[[n]][[snp_name]][["Interactions"]][["SNP"]][["Best_model"]] <- vector(mode = "list", length = snp_length)
        x <- 0
        for (snp_pair in pairs) {
          if (snp_pair[1] == snp_name) {
            x <- x + 1
            names(.master_list[[n]][[snp_name]][["Interactions"]][["SNP"]][["Best_model"]])[x] <- snp_pair[2]
            .master_list[[n]][[snp_name]][["Interactions"]][["SNP"]][["Best_model"]][[x]] <- vector(mode = "list", length = 4)
            names(.master_list[[n]][[snp_name]][["Interactions"]][["SNP"]][["Best_model"]][[x]]) <- c("Name", "Summary", "SF", "Significant")
          }
        }
      }
    }
    .master_list
  }

#' Generates a list of snp pairs to be used in further analysis steps
#'
#' Generates a list of snp pairs with (reported_model) or without (best_model)
#' the snp model.
#'
#' @param mode string, whether it should use a snp model in the interactions
#' list ("reported_model") or just the snp in question ("best_model")
#' @param snp_df dataframe, the snp with the snp interactions with the variable
#' "snp1" being the first member of the interaction and "snp2" the second. If
#' using models, "model1" should define the first member of the interaction and
#' model2 the second one.
#'
#' @return A list of character vectors, containing the character vectors the snp
#' names or snp models interactions, depending on the mode of operation selected
#' @export
generate_list_of_snp_pairs <-
  function(mode, snp_df) {

    nrows <- nrow(snp_df)
    list_of_snp_pairs <- vector(mode = "list", length = nrows)

    if (mode == "best_model") {
      for (i in seq_len(nrows)) {
        snp_1 <- snp_df$snp1[i]
        snp_2 <- snp_df$snp2[i]
        snp_pair <- c(snp_1, snp_2)
        list_of_snp_pairs[[i]] <- snp_pair
      }
    }

    if (mode == "reported_model") {
      for (i in seq_len(nrows)) {
        snp_1 <- snp_df$model1[i]
        snp_2 <- snp_df$model2[i]
        snp_pair <- c(snp_1, snp_2)
        list_of_snp_pairs[[i]] <- snp_pair
      }
    }
    list_of_snp_pairs <- unique(list_of_snp_pairs)
  }

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

# GENERIC METHODS ---------------------------------------------------------


# MISCELLANEOUS -----------------------------------------------------------
