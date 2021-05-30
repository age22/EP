
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

# Get significant results
significant_results <-
  function(list) {
  ### Generate a list of all significant results for each snp in each dataset
  for (i in seq_along(master_list)) {  # Iterating datasets
    dataset <- names(master_list)[i]
    for (n in seq_along(master_list[[i]])) { # Iterating snps
      snp <- names(master_list[[i]])[n] # Get the snp name (i.e. rs1231)
      best_model <- names(master_list[[i]][[n]][["Best_model"]]) # Get the best snp_model name (i.e. rs1231a)
      dataset_type <- ifelse(dataset == "All", "All", "Region")
      path_parts <- c(getwd(), "Output", "Main_Effects", dataset_type, dataset)
      path <- create_paths(path_parts)
      model_file <- file.path(path, paste0(snp, ".xls"))
      model <- read.xlsx(file = model_file, sheetName = best_model)
      master_list[[i]][[n]][["Main_effects"]] <- model[, c(1,4,7)] # Extracting the columns 1 (name), 3 (odds ratio) and 7 (p-value).
      if (verbose) {
        print(sprintf("Extracting significances from snp_model %s, in the %s dataset", best_model, dataset))
        cat("###############################################################################################\n")
      }
    }
  }
  master_list
}

print_significant_results <-
  function(list) {
  ma_list <- lapply(master_list, function(x) lapply(x, function(y) y[["Main_effects"]][y[["Main_effects"]][,3] < 0.05,]))
  for (i in seq_along(ma_list)) {
    print(names(ma_list)[i])
    #### Count and print total number of significant terms for each dataset
    #### (to a maximum of one term for each). This allow us to see if any
    #### variable is associated more significantly with the Diagnosis across
    #### all the dataset
    print(sort(table(unname(unlist(lapply(ma_list[[i]], function(x) x[,1])))), decreasing = T))

    #### Get the significant variables that correspond to the association between a snp_model and the diagnosis
    sign_snp <- grep("rs.*[0-9]$", unname(unlist(ma_list[[i]])), value = T)

    #### For each of these variables find the gene they belong to and if count them
    sign_genes <- vector(mode = "character", length = length(sign_snp))
    for (n in seq_along(sign_snp)) {
      snp <- sub("[a-z][0-9]$", "", sign_snp[n])
      gene <- list_of_objects[[snp]][["gene"]]
      sign_genes[n] <- gene
    }
    cat("\n")
    print(sort(table(sign_genes), decreasing = T))
    cat("\n")
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

# GENERIC METHODS ---------------------------------------------------------


# MISCELLANEOUS -----------------------------------------------------------
