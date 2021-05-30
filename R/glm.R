#' Calculates OR's, SF's and CI's from a glm.
#'
#' This function, given a logistic regression model modifies the normal output
#' of summary(model) to return also a 95% CI (upper and lower bounds) and the
#' OR/SF of the association. All credits for this code goes to Mario Cortina
#' Borja.
#'
#' @param mod1 the glm model
#' @param n.digits the number of significant digits to output
#' @param verbose whether it should print an output
#'
#' @return the summary(model) updated with the CI and the OR/SF.
#' @export
Coeff.OR2 <-
  function(mod1, n.digits=5, verbose = F) {###13/Nov/08 ###Added "verbose" argument by Armand on 2021
    ### calculates OR's, SF's and CI's from a glm.
    sum1 <- summary(mod1)
    tab <- summary( mod1)$coeff
    OR <- exp( tab[,1]); lower <- exp( tab[,1] - 1.96 * tab[,2]);
    upper <- exp(tab[,1] + 1.96*tab[,2])
    tab <- round(cbind(tab[,1], tab[,2], OR, lower, upper, tab[,4]), n.digits)
    dimnames( tab)[[2]] <- c('Estimate',
                             'Std. Error','OR','lower','upper','Pr(>|t|)')
    if (verbose) {
      print(paste('AIC',round(AIC(mod1),5), sep = ' = '), quote = F)
      print(paste(c('null deviance','null df'),
                  c(round(sum1$null.deviance, n.digits), sum1$df.null),
                  collapse = '  ', sep = '     = '),quote = F)
      print(paste(c('residual deviance','residual df'),
                  c(round(sum1$deviance,n.digits), sum1$df.residual),
                  collapse = '  ', sep = ' = '),quote = F)
    }
    tab
  }

#' Calculates AIC from a quasi-binomial
#'
#' This function inputs a quasi binomial glm model and outputs a quasi-AIC score
#' . All credits for this code goes to Mario Cortina Borja.
#'
#' @param model A generalized linear model of a quasibinomial
#'
#' @return An invisible NULL
#' @export
AIC_quasi <-
  function(model) {
    npar <- length(model$coeff); aic <- model$dev + 2*npar; print(paste('AIC quasi=',round(aic,5)),quote = F)
    invisible(NULL)
  }

# Create container for glm output
create_master_list <-
  function(datasets = DATASETS, snp_objects = list_of_objects) {

    #Creating master list structure

    # Level 1 -- length = length(DATASETS)
    master_list_names <- datasets # Defining the names of the list

    # Level 2 -- length = length(list_of_objects)
    dataset_list_names <- names(snp_objects)

    # Level 3 -- length = 4
    snp_list_names <- c("Gene", "Best_model", "Main_effects", "Interactions")

    # Level 4 -- length = 2
    interactions_list_names <- c("SNP", "Other_covariates")

    # Level 5 -- length = 3
    covariates_list_names <- c("Age75", "Sex", "E4status")

    # Level 6 -- length = 4
    variable_list_names <- c("Name", "Summary", "SF", "Significant")

    # Create recursive list
    length_vector <- c(length(DATASETS), length(list_of_objects), 4, 2, 3, 4)
    master_list <- rec_list(length_vector)

    # Name nested list and remove unnecesary levels
    names(master_list) <- master_list_names
    for (i in seq_along(master_list)) {
      names(master_list[[i]]) <- dataset_list_names
      for (j in seq_along(master_list[[i]])) {
        names(master_list[[i]][[j]]) <- snp_list_names
        master_list[[i]][[j]][["Gene"]] <- character(0)
        master_list[[i]][[j]][["Best_model"]] <- double(0)
        master_list[[i]][[j]][["Main_effects"]] <-  character(0)
        for (x in seq_along(master_list[[i]][[j]])) {
          if (!is.list(master_list[[i]][[j]][[x]])) {
            next
          }
          names(master_list[[i]][[j]][[x]]) <- interactions_list_names
          master_list[[i]][[j]][["Interactions"]][["SNP"]] <- character(0)
          for (y in seq_along(master_list[[i]][[j]][[x]])) {
            if (!is.list(master_list[[i]][[j]][[x]][[y]])) {
              next
            }
            names(master_list[[i]][[j]][[x]][[y]]) <- covariates_list_names
            for (z in seq_along(master_list[[i]][[j]][[x]][[y]])) {
              names(master_list[[i]][[j]][[x]][[y]][[z]]) <- variable_list_names
            }
          }
        }
      }
    }

    master_list # Return master_list
  }

#' Main function used to perform glm models to study association in epistasis
#' related projects
#'
#' @param .mode string, whether we want to study the "main_effects",
#' the "interaction" with covariates or the "snp" - snp interactions.
#' @param .submode if the mode chosen is the "snp" one you can further declare
#' wether you want to study the "best_model" for the best performing accoring to
#' the AIC score obtained by main_effects or just the "reported_model" if you
#' are inputting the model too.
#' @param .data  character vector, names of the datasets with all the
#' information of the inheritance models.
#' @param snps SNP_set object with all the SNPs information
#' @param covariates the covariates to control for while performing the glm.
#' @param .verbose whether it should print information to the output while
#' running.
#'
#' @return a nested list with all the information on the analyses performed
#' results
#' @export
perform_analysis <-
  function(.mode, .submode = NA, .data, snps, covariates, .verbose = verbose) {
    # Initializing progress bar
    pb <- txtProgressBar(min = 1, max = length(.data)*length(snps), style = 3)

    # Defining if we want printed output
    verbose <- .verbose

    # Creating vector to store AIC scores
    if (.mode == "main_effects") {
      master_list <- create_master_list()
      AIC_scores <- vector(mode = "double", length = 3)
      snp_model_vector <- vector(mode = "character", length = 3)
    }

    for (dataset in .data) {
      dataset_index <- which(.data == dataset) - 1
      # Counting snp number
      i <- 0
      # For each snp in the list
      if (.mode == "snp") {
        names(snps) <- sapply(snps, function(x) paste(x[1], x[2], sep = "_"))
      } else if (.mode == "snp_model") {
        names(snps) <- sapply(snps, function(x) paste(x[1], x[2], sep = "_"))
      }
      for (snp in names(snps)) {
        if (dataset == "America") {
          if (snp %in% not_in_america | snp == "rs11754661") {
            next
          }
        }
        if (.mode == "snp") {
          if (.submode == "best_model") {
            covariate <- sub(pattern = "^rs[0-9]+_", replacement = "", x = snp)
            snp <- sub(pattern = "_rs[0-9]+$", replacement = "", x = snp)
          } else if (.submode == "reported_model") {
            covariate_model <- sub(pattern = "^rs[0-9]+[a-z]_", replacement = "", x = snp)
            snp_model <- sub(pattern = "_rs[0-9]+[a-z]$", replacement = "", x = snp)
            snp <- sub(pattern = "[a-z]$", replacement = "", x = snp_model)
            covariate <- sub(pattern = "[a-z]$", replacement = "", x = covariate_model)
          }
        }
        # Find in which column of the dataset they appear
        snp_index <- which(colnames(get(dataset)) == snp)
        # And keep count of the number of snps we have done yet to present it in a progress bar
        i <- i + 1
        setTxtProgressBar(pb, i + dataset_index*length(snps))
        if (.mode == "main_effects") {
          n <- 1
          while (n < 4) {
            snp_model <- colnames(get(dataset))[snp_index + n]
            result_ME <- glm_loop(mode = .mode, data = dataset, .snp_model = snp_model, .snp = snp, number = i, .covariates = covariates)
            AIC_scores[n] <- result_ME[[1]]
            snp_model_vector[n] <- snp_model
            n <- n + 1
          }
          #Adding gene information
          gene <- snps[[snp]][["gene"]]
          master_list[[dataset]][[snp]][["Gene"]] <- gene
          # Adding best AIC score (minimum one)
          master_list[[dataset]][[snp]][["Best_model"]] <- min(AIC_scores)
          # Finding out to which inheritance model does that AIC score belong
          index_AIC <- which(AIC_scores == min(AIC_scores))
          names(master_list[[dataset]][[snp]][["Best_model"]]) <- snp_model_vector[index_AIC]
        }

        if (.mode == "interaction") {
          for (covariate in covariates) {
            if (covariate == "Centre") {
              next
            }
            snp_model <- names(master_list[[dataset]][[snp]][["Best_model"]])
            results_INT <- glm_loop(mode = "interaction", data = dataset, number = i, .covariate = covariate, .covariates = covariates, .snp_model = snp_model, .snp = snp)
            master_list[[dataset]][[snp]][["Interactions"]][["Other_covariates"]][[covariate]][["Name"]] <- results_INT[[1]]
            master_list[[dataset]][[snp]][["Interactions"]][["Other_covariates"]][[covariate]][["Summary"]] <- results_INT[[2]]
            master_list[[dataset]][[snp]][["Interactions"]][["Other_covariates"]][[covariate]][["SF"]] <- results_INT[[3]]
            master_list[[dataset]][[snp]][["Interactions"]][["Other_covariates"]][[covariate]][["Significant"]] <- results_INT[[4]]
          }
        }

        if (.mode == "snp") {
          if (.submode == "best_model") {
            snp_model <- names(master_list[[dataset]][[snp]][["Best_model"]])
            covariate_model <- names(master_list[[dataset]][[covariate]][["Best_model"]])
            results_SNP <- glm_loop(mode = "snp", submode = "best_model", data = dataset, .covariate = covariate_model, .covariates = covariates, .snp_model = snp_model, .snp = snp, number = i)
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Best_model"]][[covariate]][["Name"]] <- results_SNP[[1]]
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Best_model"]][[covariate]][["Summary"]] <- results_SNP[[2]]
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Best_model"]][[covariate]][["SF"]] <- results_SNP[[3]]
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Best_model"]][[covariate]][["Significant"]] <- results_SNP[[4]]

          }
          if (.submode == "reported_model") {
            results_SNP <- glm_loop(mode = "snp", submode = "reported_model", data = dataset, .covariate = covariate_model, .covariates = covariates, .snp_model = snp_model, .snp = snp, number = i)
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Reported"]][[snp_model]][[covariate]][[covariate_model]][["Name"]] <- results_SNP[[1]]
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Reported"]][[snp_model]][[covariate]][[covariate_model]][["Summary"]] <- results_SNP[[2]]
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Reported"]][[snp_model]][[covariate]][[covariate_model]][["SF"]] <- results_SNP[[3]]
            master_list[[dataset]][[snp]][["Interactions"]][["SNP"]][["Reported"]][[snp_model]][[covariate]][[covariate_model]][["Significant"]] <- results_SNP[[4]]
          }
        }
      }
    }
    close(pb)
    master_list
  }

glm_loop <-
  function(mode = .mode, submode = NA, data = .data, .covariate = covariate, .covariates = covariates, .snp_model = snp_model, .snp = snp, number, .model = model, .verbose = verbose) {

    # Print relevant information
    if (.verbose) {
      cat("\n####################################################\n")
      print(sprintf("Dataset -> %s", data))
      model <- substr(.snp_model, nchar(.snp_model), nchar(.snp_model))
      print(sprintf("snp number %s, model %s", number, model))
      print(sprintf("snp_model -> %s", .snp_model))
      if (mode == "snp") {
        model <- substr(.covariate, nchar(.covariate), nchar(.covariate))
        print(sprintf("snp number %s, model %s", number, model))
        print(sprintf("snp_model -> %s", .covariate))
      }
      if (mode == "interaction") {
        print(.covariate)
      }
    }

    # Creating formula for glm
    if (mode == "interaction") {
      .covariates <- setdiff(.covariates, .covariate)
    }

    snp_term <- ifelse(mode == "interaction" | mode == "snp", paste0(.snp_model, "*", .covariate), .snp_model)
    predictors <- c(.covariates, snp_term)

    if (data == "America") {
      predictors <- setdiff(predictors, "Centre")
    }

    formula <- as.formula(paste("Diag", paste(predictors, collapse = " + "), sep = " ~ "))

    # Defining arguments to pass to the glm function
    args <-  list(formula = formula, family = quasibinomial("logit"), data = as.name(data))

    # Execute glm function
    linear_model <- do.call(glm, args)
    OR_summary <- Coeff.OR2(linear_model)


    if (mode == "main_effects") {

      # Calculate AIC score for the model
      AIC <- capture.output(AIC_quasi(linear_model))
      result_ME <- list(AIC)

      # Generate path where the excel file will be stored
      dataset_type <- ifelse(data == "All", "All", "Region")
      path_parts <- c(output, "Main_Effects", dataset_type, data)
      path <- create_paths(path_parts)

      # Write to excel
      xlsx::write.xlsx2(AIC, file = file.path(path, paste0(.snp, ".xls")), sheetName = paste0("AIC", "_", .snp_model), append = T, sep = "\t")
      xlsx::write.xlsx2(OR_summary, file = file.path(path, paste0(.snp, ".xls")), sheetName = .snp_model, append = T, sep = "\t")
    }

    if (mode == "interaction" | mode == "snp") {

      # Extracting interaction term from the summary
      if (mode == "interaction") {
        pattern_12 <- paste0(.snp_model, '[0-9]:', .covariate, '.*$')
      }

      if (mode == "snp") {
        pattern_12 <- paste0(.snp_model, '[0-9]:', .covariate, '[0-9]$')
      }

      interaction_indexes <- grep(pattern_12, rownames(OR_summary))
      interaction_terms <- grep(pattern_12, rownames(OR_summary), value = T)

      terms <- vector(mode = "character", length = length(interaction_terms))
      SF <- vector(mode = "double", length = length(interaction_terms))
      is_significant <- vector(mode = "logical", length =  length(interaction_terms))

      for (level in levels(get(data)[[.covariate]])) {
        my_table <- table(get(data)[[.snp_model]][get(data)[[.covariate]] == level])
        for (count in my_table) {
          if (count < 2) {
            if (.verbose) {
              print(sprintf("Variable %s for snp %s doesn't have enough counts for level %s to perform an analysis", .covariate, .snp_model, level))
              print(my_table)
            }
            assign(ifelse(mode == "interaction", "result_INT", "result_SNP"), list("Not enough counts", list(NULL), 0, FALSE))
          }
        }
      }

      for (i in seq_along(interaction_terms)) {
        interaction_index <- interaction_indexes[i]
        interaction_term <- interaction_terms[i]

        # Extracting covariate term from the summary
        pattern_1 <- sub("^.*:", "", interaction_term)

        # Extracting snp_model term from the summary
        pattern_2 <- sub(":.*$", "", interaction_term)

        # Extracting SF from the summary
        SF[i] <- OR_summary[interaction_index, 3]
        terms[i] <- interaction_term
        is_significant[i] <- ifelse(OR_summary[interaction_index, 6] < 0.05, TRUE, FALSE)
        assign(x = ifelse(mode == "interaction", "result_INT", "result_SNP"), value = list(terms, OR_summary, SF, is_significant))

        # Generate path where the excel file will be stored
        dataset_type <- ifelse(data == "All", "All", "Region")
        if (mode == "interaction") {
          path_parts <- c(output, "Interaction", dataset_type, data)
        } else {
          path_parts <- c(output, "SNP-SNP")
          if (submode == "reported_model") {
            path_parts <- c(path_parts, "Reported", dataset_type, data)
          } else if (submode == "best_model") {
            path_parts <- c(path_parts, "Best_model", dataset_type, data)
          }
        }
        path <- create_paths(path_parts)

        # Write to excel
        write.xlsx2(OR_summary, file = file.path(path, paste0(.snp_model, ".xls")), sheetName = paste0("OR", "_", pattern_1, "_", pattern_2), append = T, sep = "\t")
        write.xlsx2(is_significant, file = file.path(path, paste0(.snp_model, ".xls")), sheetName = paste0("SF", "_", pattern_1, "_", pattern_2), append = T, sep = "\t")
      }
    }

    # Return value to main function
    if (mode == "main_effects") {
      result <- result_ME
    } else if (mode == "interaction") {
      result <- result_INT
    } else {
      result <- result_SNP
    }
    result
  }
