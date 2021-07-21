#' Prints main effects results
#'
#' This function prints main effects results, applying multiple testing
#' correction or not. In particular prints the snps with their p-values and the
#' genes to which they correspond.
#'
#' @param list list, containing the results in the analysis (i.e. master_list)
#' @param corrected boolean, whether it should correct for multiple testing.
#' Default is false.
#'
#' @return NUll (invisible)
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
      p_values <- stats::p.adjust(p_values, method = "bonferroni")
      significant <- which(p_values < 0.05)
      p_values <- p_values[significant]
      snps <- snps[significant]
      is_best_model <- sapply(snps, function(snp_model_lvl) {
                              snp_model <- sub("[0-2]$", "", snp_model_lvl);
                              snp <- sub("[a-z]$", "", snp_model);
                              best_model <- list$All[[snp]]$Best_model;
                              best_model <- names(best_model);
                              identical(snp_model, best_model)
      })
      if (length(is_best_model) > 0) {
        p_values <- p_values[is_best_model]
        snps <- snps[is_best_model]
      }
    }

    names(p_values) <- snps

    if (length(p_values > 0)) {
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

  } else {
    cat("No significant interactions found")
  }

    invisible(NULL)
}

#' Print interaction between SNP and covariates Age and Sex.
#'
#' This function prints covariate interaction results, applying for multiple
#' testing correction or not. In particular prints the snp-covariate interaction
#' with its p-value and the genes to which they correspond.
#'
#' @param list list, containing the results in the analysis (i.e. master_list)
#' @param covariates vector, with the covariates that we want to study
#' @param corrected boolean, whether it should correct for multiple testing.
#' Default is false.
#'
#' @return NULL (invisible)
#' @export
print_INT <-
  function(list, covariates, corrected = F) {

    # Getting just the significant Interaction (INT) p-value column
    if (!corrected) {

    sig_INT <- lapply(list$All, function(snp_list) {
      interactions <- snp_list$Interactions$Other_covariates[-c(1,3)]
      lapply(interactions, function(int_list){
        interaction <- int_list$Summary;
        interaction[interaction[,6] < 0.05,][,6]})})

    # Remove snp info and just keep the pvalues for each significant term
    terms <- lapply(unname(sig_INT), function(x) {unname(x)})

    # Getting all Interaction (INT) p-value column
    } else {

      sig_INT <- lapply(list$All, function(snp_list) {
        interactions <- snp_list$Interactions$Other_covariates[-c(1,3)]
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
      p_values <- stats::p.adjust(p_values, method = "bonferroni")
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
        browser()
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

#' Print SNP-SNP interactions results.
#'
#' This function prints SNP-SNP interaction results, applying for multiple
#' testing correction or not. In particular prints the snp-snp interaction
#' with its p-value and the genes to which they correspond.
#'
#' @param master_list list, containing the results in the analysis
#' @param list_of_snp_pairs list, containing the SNP-SNP interactions
#' to look for
#' @param list_of_APOE list, containing the SNP-APOE4 interactions to look for
#' @param corrected boolean, whether it should correct for multiple testing.
#' Default is false.
#'
#' @return NULL (invisible)
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
      p_values <- stats::p.adjust(p_values, method = "bonferroni")
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
      cat("No significant interactions found\n")
    }
    invisible(NULL)
  }
