
# PREPARATION -------------------------------------------------------------


# Constructor for the "SNP" object
new_SNP <-
  function(gene, id, geno_count, geno_freq, allele_count, allele_freq, minor_allele,
  major_allele, imput_geno_count, imput_geno_freq, imput_allele_count, imput_allele_freq,
  imput_minor_allele, imput_major_allele) {

  stopifnot(
    is.character(gene), is.character(id), is.integer(geno_count),
    is.double(geno_freq), is.double(allele_count), is.double(allele_freq),
    is.character(minor_allele), is.character(major_allele),
    is.integer(imput_geno_count), is.double(imput_geno_freq),
    is.double(imput_allele_count), is.double(imput_allele_freq),
    is.character(imput_minor_allele), is.character(imput_major_allele))

  structure(
    list(
      gene = gene,
      id = id,
      geno_count = geno_count,
      geno_freq = geno_freq,
      allele_count = allele_count,
      allele_freq = allele_freq,
      minor_allele = minor_allele,
      major_allele = major_allele,
      imput_geno_count = imput_geno_count,
      imput_geno_freq = imput_geno_freq,
      imput_allele_count = imput_allele_count,
      imput_allele_freq = imput_allele_freq,
      imput_minor_allele = imput_minor_allele,
      imput_major_allele = imput_major_allele
    ),
  class = "SNP"
  )
}

# List of objects creator function
generate_object <-
  function(exists, snps, ...) {
  # Checking that the exists arguments is boolean (Either TRUE or FALSE)
  if (!is.logical(exists)) {
    stop("The 'exists' argument must be a boolean")
  }

  # If a list with all SNP objects doesn't exist (= FALSE), create an empty list
  if (!exists) {
    list_of_objects <- vector(mode = "list", length = length(snps))
  }

  # For each snp
  for (i in seq_along(snps)) {
    # Get snp name
    snp <- snps[i]

    # Generate SNP object attributes and output them in a list (SNP_attributes)
    SNP_attributes <- create_SNP_attributes(snp = snp, ...)

    # Then assign each attribute to its variable
    geno_count <- SNP_attributes[[1]]
    geno_freq <- SNP_attributes[[2]]
    allele_count <- SNP_attributes[[3]]
    allele_freq <- SNP_attributes[[4]]
    minor_allele <- SNP_attributes[[5]]
    major_allele <- SNP_attributes[[6]]

    # If a list with all SNP objects doesn't exist (= FALSE), create a SNP
    # object for the current SNP
    if (!exists) {
      # Split the gene_snp string into gene and snp
      gene_snp_pair <- strsplit2vector(snp, pattern = "_")

      # Creating a SNP object and defining the gene and id attributes
      object <- new_SNP(gene = gene_snp_pair[1],
                        id = gene_snp_pair[2],

                        geno_count = geno_count,
                        imput_geno_count = integer(),

                        geno_freq = geno_freq,
                        imput_geno_freq = double(),

                        allele_count = allele_count,
                        imput_allele_count = double(),

                        allele_freq = allele_freq,
                        imput_allele_freq = double(),

                        minor_allele = minor_allele,
                        imput_minor_allele = character(),

                        major_allele = major_allele,
                        imput_major_allele = character()
                        )

    # The list_of_objects already exists
    } else {

      # As the object has already been created, modify existing SNP object
      object <- list_of_objects[[snp]]
      object[["imput_geno_count"]] <- geno_count
      object[["imput_geno_freq"]] <- geno_freq
      object[["imput_allele_count"]] <- allele_count
      object[["imput_allele_freq"]] <- allele_freq
      object[["imput_minor_allele"]] <- minor_allele
      object[["imput_major_allele"]] <- major_allele

    }

    list_of_objects[[i]] <- object # Adding object created/modified to the
                                   # list of objects

    if (!exists) {
      # Defining object name as its snp id
      names(list_of_objects)[i] <- object$id
    }

  }
  # Adding personalized class to the list so methods print.SNP_set and
  # str.SNP_set for their respective generics can be defined.
  class(list_of_objects) <- 'SNP_set'

  # Return filled list_of_objects
  list_of_objects
}

create_SNP_attributes <-
  function(snp, dataset) {
    ## Generating genotypic counts by doing a summary table of the
    ## genotypes in the control cases for a specific snp
    controls_subset <- dataset[dataset$Diag == "Control",]
    snp_controls <- controls_subset[[snp]]
    geno_count <- table(snp_controls)

    ## Calculate genotype frequencies
    geno_freq <- 100*geno_count/sum(geno_count)
    geno_freq <- round(geno_freq, digits = 3)
    ## Calculate allele counts
    allele_count <- calculate_allele_counts(geno_count, dataset, snp)

    ## Calculate allele frequencies
    allele_freq <- 100*allele_count/sum(allele_count)
    allele_freq <- round(allele_freq, digits = 3)

    ## Calculate minor and major allele.
    minor_allele <- names(which.min(allele_count))
    major_allele <- names(which.max(allele_count))

    ## Returning list with SNP_attributes
    list(geno_count,
         geno_freq,
         allele_count,
         allele_freq,
         minor_allele,
         major_allele)
  }

calculate_allele_counts <-
  function(genotype_count, dataset, snp) {
    ### Allocating container variables that will contain the allele names
    ### obtained from the genotype names as well as their respective counts.
    ### Each genotype will result in 2 alleles. Then, given 3 genotypes, we will
    ### obtain 6 alleles (allowing repetition) with their respective counts.
    allele_names <- character(6)
    allele_values <- numeric(6)

    ### If the controls for a certain snp don't show all 3 genotypes modify
    ### the variable 'how_many_genotypes' as we are expecting
    ### (2 * how_many_genotypes) alleles instead.
    how_many_genotypes <- length(genotype_count)
    if (how_many_genotypes != 3) {
      allele_names <- character(2*how_many_genotypes)
      allele_values <- numeric(2*how_many_genotypes)

      # (Optional) Print information about which SNP doesn't have all 3
      # genotypes for the controls
      if (verbose) {
        # Number of genotypes for control + alzheimer cases
        how_many_total_genotypes <- length(table(dataset[[snp]]))
        info <- c("For the control subset snp %s there are only %s genotypes",
                  "available instead of the %s found in the total dataset",
                  "(including diagnosed cases). They are the following ones:")
        info <- paste(info, collapse = " ")
        print(sprintf(info, snp, how_many_genotypes, how_many_total_genotypes))
        print(genotype_count)
      }
    }

    ### Generate a vector with the allele names (i.e "A", "T", "G" or "C")
    ### in a snp and another vector with the counts (i.e 50, 25, 75...).
    for (x in seq_along(genotype_count)) {
      genotype <- names(genotype_count)[x]
      count <- genotype_count[x]

      #### Split the genotype into each of its alleles and transform the list
      #### obtained into a character vector
      alleles <- strsplit2vector(genotype, pattern = "")

      # Each genotype position (x) results in two allele positions,
      # position 2x-1 and position 2x in the allele_names list.
      # Example: genotype position 2 results in allele
      # position 3 (2*2-1) and allele position 4 (2*2).
      y <- 2*x - 1

      for (allele in alleles) {
        allele_names[y] <- allele
        allele_values[y] <- count
        y <- 2*x
      }
    }

    ### Merge the allele_counts vector with the allele_names one, split the
    ### allele_counts according to allele name and then sum the counts inside of
    ### each allele. i.e A = 25, A = 25, A = 50, T = 50, T = 75, T = 75 gets
    ### split as A = c(25, 25, 50) and T = c(50, 75, 75) and then summed as
    ### as A = 100 and T = 200.
    names(allele_values) <- allele_names
    allele_count <- split(allele_values, allele_names)
    allele_count <- sapply(allele_count, sum)
  }

genotype_imputated_df <-
  function(list_of_objects, df, thousand_genomes, verbose, ...) {

    # Container of snp names with non-matching reference and genotyped alleles
    wrong <- vector(mode = "list", length = length(list_of_objects))

    for (i in seq_along(list_of_objects)) {
      object <- list_of_objects[[i]]

      snp <- object$id
      major_allele <- object$major_allele
      minor_allele <- object$minor_allele
      genotyped_alleles <- c(major_allele, minor_allele)

      ### Calculating the mean of all the individuals for the allele dosage at a certain snp.
      snp_mean <- mean(df[[snp]], na.rm = TRUE)

      if (verbose) {
        cat("------------------------------------------------------------------------------------\n")
        print(sprintf("The mean of all rows from the snp %s is %s", snp, round(snp_mean, digits = 3)))
        print(sprintf("Confirming that there is agreement in the reference chromosomal strand between the imputed data and the genotyped data for the snp %s.", snp))
      }

      #### If the reference or alternative allele are the same as the major or minor allele it
      #### probably means it is in strand agreement, else it may be the complementary strand
      snp_id <- thousand_genomes$snp_id == snp
      ref_allele <- thousand_genomes[snp_id,]$reference
      alt_allele <- thousand_genomes[snp_id,]$alternative
      `1000_geno_alleles` <- c(ref_allele, alt_allele)

      #### If none of the strands, neither the forward nor the reverse strand
      #### agrees with the assigned genotype alleles, it will check
      #### which alleles show such disagreement in order to look at it further
      set1 <- genotyped_alleles
      set2 <- `1000_geno_alleles`
      set3 <- complement_alleles(set2)
      if (setequal(set1,set2) == FALSE & setequal(set1, set3) == FALSE) {
        wrong[[i]] <- snp

        if (verbose) {
          print(sprintf("The chromosomal reference alleles and the genotyped data alleles don't match for snp %s", snp))
          print(set1)
          print(set2)
          print(set3)
        }
        next
      }
      ### Allele dosage data is recoded into one of the three genotypes possible
      ### in which a value higher than 1.5 means that the individual is homozygote
      ### for allele 1, if lower than 0.5 homozygote for allele 2 and if in between
      ### heterozygote


      allele_1 <- ifelse(ref_allele %in% genotyped_alleles,
                         ref_allele,
                         complement_alleles(ref_allele))

      allele_2 <- ifelse(alt_allele %in% genotyped_alleles,
                         alt_allele,
                         complement_alleles(alt_allele))

      if (exists("weird")) {
        if (snp %in% weird) {
          a1 <- allele_1
          a2 <- allele_2
          allele_1 <- a2
          allele_2 <- a1
        }
      }

      if (verbose) {
        print(allele_1)
        print(allele_2)
      }

      homozygote_1 <- paste0(allele_1, allele_1)
      homozygote_2 <- paste0(allele_2, allele_2)
      heterozygote <- paste0(allele_1, allele_2)

      df[[snp]][df[[snp]] > 1.5] <- homozygote_1
      df[[snp]][df[[snp]] < 0.5] <- homozygote_2
      df[[snp]][df[[snp]] >= 0.5 & df[[snp]] <= 1.5] <- heterozygote
    }

    wrong <- unlist(wrong)
    wrong <- wrong[is.null(wrong) == FALSE]
    if (verbose) {
      cat("\n")
      if (length(wrong) != 0) {

        print(sprintf("%s snp not found in the thousand genomes dataset", wrong))

      } else {

        print("All snps were found in the thousand genomes dataset")

      }
    }
    df
  }

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


# Set up log reg analysis
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

# Calculates AIC from a quasi-binomial
AIC.quasi <-
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

# Perform GLM model
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
            snp_model <- colnames(get(.data))[snp_index + n]
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
      AIC <- capture.output(AIC.quasi(linear_model))
      result_ME <- list(AIC)

      # Generate path where the excel file will be stored
      dataset_type <- ifelse(data == "All", "All", "Region")
      path_parts <- c(output, "Main_Effects", dataset_type, data)
      path <- create_paths(path_parts)

      # Write to excel
      write.xlsx2(AIC, file = file.path(path, paste0(.snp, ".xls")), sheetName = paste0("AIC", "_", .snp_model), append = T, sep = "\t")
      write.xlsx2(OR_summary, file = file.path(path, paste0(.snp, ".xls")), sheetName = .snp_model, append = T, sep = "\t")
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

print.SNP <-
  function(SNP) {
    cat("######", SNP$id, "######")
    for (i in seq_along(SNP)) {
      attribute <- names(SNP)[i]
      value <- SNP[[i]]
      name <- names(SNP[[i]])
      exceptions <- c("gene", "id", "minor_allele", "major_allele", "imput_minor_allele", "imput_major_allele")
      if (length(value) != 0) {
        if (attribute %in% exceptions) {
          cat(sprintf("\n\n[%s] %s", attribute, value))
        } else {
          cat(sprintf("\n\n[%s]", attribute))
          cat(sprintf(" %s = %s;", name, value))
        }
      }
    }
  }

# print() and str() generics methods for SNP_set objects (i.e. list_of_objects)
print.SNP_set <-
  function(SNP_set_object) {
    cat("This is just a convenience function, invocated when print() is used, to provide an overview of the genes and snps present in this list. To view its contents, use list_of_objects$your_snp_id\n\n\n")
    snps <- names(SNP_set_object)
    genes <- unique(unname(sapply(SNP_set_object, function(x) x$gene)))
    print(list('snps' = snps, 'genes' = genes))
    cat(sprintf("\nThe list has a total of %s snps and %s genes", length(snps), length(genes)))
  }

str.SNP_set <-
  function(SNP_set_object) {
    cat("This is just a convenience function, invocated when str() is used, to provide an overview of the list elements. To view its contents, use list_of_objects$your_snp_id\n\n")
    cat("An example of the structure of the objects in the list is the following one:\n\n")
    print(ls.str(list_of_objects[1]))
  }


# MISCELLANEOUS -----------------------------------------------------------

load_packages <-
  function(packages) {
    invisible(suppressPackageStartupMessages(
        lapply(packages, library, character.only = TRUE)))
  }

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
  }

complement_alleles <-
  function(sequence) {

  # Creating complement function
  complement_function <- function(allele) {
    if (allele == "G") {
      complement <- "C"
    }
    if (allele == "C") {
      complement <- "G"
    }
    if (allele == "A") {
      complement <- "T"
    }
    if (allele == "T") {
      complement <- "A"
    }
    complement
  }
  # Using complement function on sequence elements
  complementary <- sapply(sequence, complement_function)
  complementary <- unname(complementary)

  # Returning complementary sequence obtained
  complementary
}

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

rec_list <-
  function(x) {
  if (length(x) == 1) {
    vector("list", x)
  } else {
    lapply(1:x[1], function(...) rec_list(x[-1]))
  }
  }

strsplit2vector <-
  function(string, pattern) {
    unlist(strsplit(x = string, split = pattern))
  }
