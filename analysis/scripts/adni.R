#//////////////////////////////PREPARATION//////////////////////////////////////

# 1. Genotyped data set preparation ------------------------------------------

  ## Clearing global environment
  rm(list = ls())

  ## Defining input, output, scripts and functions path shortcuts.
  input <- file.path(getwd(), "analysis/data/Input")
  output <- file.path(getwd(), "analysis/data/Output")
  scripts <- file.path(getwd(), "analysis/scripts")
  functions <- file.path(scripts, "functions.R")

  ## Deleting the output folder and all the files and folders it contains
  ## from the previous run
  unlink(output, recursive = TRUE)

  ## Importing functions from the "functions" script
  source(functions)

  ## Loading packages

    ### Naming packages names
    my_packages <- c("psych", "stringr", "xlsx", "ggplot2")

    ### Loading named packages silently
    load_packages(my_packages)

  ## Check if it is run as a script (i.e. from command line) and if so, get
  ## valid arguments inputted by the user and assigning them as TRUE while
  ## keeping the rest as FALSE

    ### THESE ARE THE ARGUMENTS THAT THESE R SCRIPT CAN TAKE FROM THE USER
    ### (DEFAULT IS TRUE)
    ### --verbose means whether to output to the terminal
    ### --imputated is if the program should expect an imputated dataset
    ### --EP performs needed data correction before analysis in the Epistasis
    ###   Project dataset

    if (!interactive()) {

      ### Get arguments inputted by the user
      arguments <- commandArgs(trailingOnly = TRUE)

      ### Vector of allowed options
      options <- c("--verbose", "--imputated", "--EP")

      ### Assigning the valid arguments inputted by the user as TRUE and the
      ### rest as FALSE.
      assign_arguments(arguments, options)

    } else {

      # If you want to modify the behaviour of the script while running it
      # interactively, change the boolean value of these options
      verbose <- TRUE
      imputated <- TRUE
      EP <- TRUE
      ADNI <- FALSE

    }

  ## Reading and storing genotyped input data

    ### Printing to output
    cat("\nPerforming genotyped dataset preparation...\n\n\n")

    ### Reading csv file (";" as separator) and storing it as a dataframe
    encoded_NA <- c("00", "", "???", "-9")
    df <- read.csv2(file.path(input, "genotyped.csv"),
                    na.strings = encoded_NA)

    ### Fixing several issues from the snps in the Epistasis project dataset
    if (EP) {
      df$DHFR_rs70991108_INDEL <- NULL # Deleting indel from the dataset
      df$PPARG_rs709149 <- NULL
      df$MS4A4E_rs670139 <- NULL
      # swap1 <- df$PPARG_rs709149 # Fixing suspected swap
      # swap2 <- df$MS4A4E_rs670139
      # df$PPARG_rs709149 <- swap2
      # df$MS4A4E_rs670139 <- swap1
      # rm(list = c("swap1", "swap2"))
    }

  ## Subsetting and extracting input data into useful data structures

    ### Reading the snps to subset from snps.txt file and defining the vector
    ### type as character
    vector_of_gene_snps <- scan(file.path(input, "snps.txt"),
                                what = "character",
                                quiet = !verbose)

    ### Choosing variables of interest
    variables <- c("ID", "Diag", "SexLett", "Age", "Centre", "LGC_E4.")

    ### Subsetting the dataframe with the desired variables and SNPs
    df <- df[, c(variables, vector_of_gene_snps)]

    if (ADNI) {
      adni <- read.csv(file.path(input, "adni.csv"))
      adni_snps <- colnames(adni)[7:ncol(adni)]
      variables <- c("ID", "Diag", "SexLett", "Age", "Centre", "LGC_E4.")
      adni <- adni[,c(variables, adni_snps)]
      not_in_america <- sub(setdiff(vector_of_gene_snps, adni_snps), pattern = ".*_", replacement = "")
      df <- merge(df, adni, all = T)
      df <- df[, c(variables, vector_of_gene_snps)]
    }

    ### Creating a list of objects that for each SNP contains information about
    ### its name, gene it belongs, genotype counts, genotype frequencies,
    ### allele frequencies and counts, minor allele and major allele etc.
    list_of_objects <- generate_object(exists = exists('list_of_objects'),
                                       dataset = df,
                                       snps = vector_of_gene_snps)

    ### Creating a vector of names of the snps available in the dataset
    ### (same as vector_of_gene_snps but without the gene in the name)
    vector_of_snps <- sapply(list_of_objects, function(x) x$id)
    vector_of_snps <- unname(vector_of_snps)

    ### Converting categorical variables from character to factor data type.
    df[-c(1,4)] <- lapply(df[-c(1,4)], as.factor)

  ## Remove used variables
  objects_to_be_removed <- c("variables", "my_packages", "vector_of_gene_snps",
                             "scripts", "functions")

  rm(list = objects_to_be_removed)


# 2. Imputed dataset preparation ----------------------------------------------

  if (imputated) {
    ## Printing to output
    cat("Performing imputated dataset preparation...\n\n\n")

    ## Reading and storing imputated input data
    encoded_NA <- c("#NULL!", "NA")
    imput_df <- read.csv2(file.path(input, "imputated.csv"),
                          na.strings = encoded_NA)

    ## Subsetting and reformatting input data
    imput_variables <- c("id", "DiagName", "sex", "Age.to.Use",
                         "Age75", "E4status")
    imput_df <- imput_df[,c(imput_variables, vector_of_snps)]

    ## Importing data from the reference panel of the 1000 genomes project
    ## phase 1 version 3 in order to be able to recode the allele dosage data
    ## into genotypes
    thousand_genomes <- read.table(file.path(input, "snps_1000genomes.txt"))
    colnames(thousand_genomes) <- c("chromosome", "location", "snp_id",
                                    "reference", "alternative")

    if (EP) {
      ## Renaming thousand genomes snp_id for alternative id names which are
      ## the ones used in the imputated dataset
      select_snp <- thousand_genomes$snp_id == "rs116803374"
      thousand_genomes[select_snp,]$snp_id <- "rs242557"

      select_snp <- thousand_genomes$snp_id == "rs115050875"
      thousand_genomes[select_snp,]$snp_id <- "rs2471738"

      ## snp rs1052533 not found in thousand_genomes. Directly genotyped in
      ## imput_df maybe? So we add the reference genome information from dbSNP
      ## to be able to convert the "imputation" made to genotype data.

      index <- nrow(thousand_genomes) + 1 # We'll add a new snp to the thousand
                                          # genomes dataset extracted from
                                          # current (1st quarter of 2021) data
                                          # of dbSNP

      thousand_genomes[index, "chromosome"] <- 19
      thousand_genomes[index, "location"] <- 50378568
      thousand_genomes[index, "snp_id"] <- "rs1052533"
      thousand_genomes[index, "reference"] <- "G"
      thousand_genomes[index, "alternative"] <- "A"
    }

    ## Converting imputation data type to numeric.

    imput_df[, 7:length(imput_df)] <- lapply(imput_df[,7:length(imput_df)],
                                            function(x) as.numeric(x))

    if (EP) {
      ## SNPs that seemed to have reversed genotypic frequencies
      weird <- c("rs5749131", "rs1801198", "rs7289553", "rs5997711", "rs326120",
                 "rs162040", "rs2228604", "rs10745354", "rs669340", "rs1611115",
                 "rs242557", "rs2471738", "rs4800488", "rs2695121", "rs1187323",
                 "rs1545285", "rs2069442", "rs3931914", "rs1800896", "rs670139",
                 "rs709149", "rs13022344")

      ## Recoding SNP continuous variable into genotypes in the imputated dataframe.
      imput_df <- genotype_imputated_df(list_of_objects = list_of_objects,
                                        df = imput_df,
                                        thousand_genomes = thousand_genomes,
                                        verbose = verbose,
                                        weird = weird)

    } else {

      ## Recoding SNP continuous variable into genotypes in the imputated dataframe.
      imput_df <- genotype_imputated_df(list_of_objects = list_of_objects,
                                        df = imput_df,
                                        thousand_genomes = thousand_genomes,
                                        verbose = verbose)
    }
    ## Remove used variables
    objects_to_be_removed <- c("index", "weird", "imput_variables", "encoded_NA")
    suppressWarnings(rm(list = objects_to_be_removed))
  }




# 3.	Merge imputed and genotyped data sets--------

  ## Uniformizing data and reformatting data and preliminary QC steps for the genotyped dataset.

    ### Renaming variable names in the genotyped dataset.
    colnames(df) <- c("ID", "Diag", "Sex", "Age_to_use", "Centre", "E4status", vector_of_snps)

    ### As we are studying LOAD, subsetting individuals higher or equal than 60 years old in the genotyped dataset
    df <- df[df$Age_to_use >= 60,]

    ### Merging variables

      #### E4status

        ##### Renaming the levels of the E4status binary values
        levels(df$E4status) <- c("E4-", "E4+")

      #### Age75

        ##### Adding the Age75 variable to the genotyped dataset and renaming the
        ##### levels of other categorical variables in the datasets.
        df$Age75 <- ifelse(df$Age_to_use >= 75, ">75", "<75")


  ## Uniformizing data and reformatting data and preliminary QC steps for the imputated dataset.

  if (imputated) {

    cat("Uniformizing and merging the datasets...\n\n\n")

    ### Converting data type from characters to factors
    imput_df[-c(1,4)] <- lapply(imput_df[-c(1,4)], as.factor)

    ### Renaming variable names in the imputated dataset in order to simplify merging later
    colnames(imput_df) <- c("ID", "Diag", "Sex", "Age_to_use", "Age75", "E4status", vector_of_snps)

    ### As we are studying LOAD, subsetting individuals higher or equal than 60 years old in the imputated dataset
    imput_df <- imput_df[imput_df$Age_to_use >= 60,]

    ### Merging variables

      #### Diagnosis

        ##### Renaming the levels of the Diagnostic categorical variable for the imputated dataset
        levels(imput_df$Diag) <- c("Control", "AD")

      #### Sex

        ##### Eliminating some individuals of unclear sex as there is a factor termed "1" in
        ##### the otherwise Male or Female levels of the sex variable and renaming the
        ##### levels of the sex variable for the <- rotterdam dataset
        imput_df <- imput_df[imput_df$Sex != 1,]
        imput_df$Sex <- factor(imput_df$Sex) # This is done to reset the levels to two instead of the previous (wrongly assigned) three
        levels(imput_df$Sex) <- c("Male", "Female")

      #### Age75

        levels(imput_df$Age75) <- c("<75", ">75")

      ####  Centre

        ##### Assigning a new variable to the Rotterdam Dataset to trace the
        ##### Rotterdam individuals as it was done with the other centres in the
        ##### genotyped database
        imput_df$Centre <- "ROTTERDAM"

  ## Adding the attributes genotype counts, genotype frequencies, allele counts
  ## allele frequency, minor allele and major allele to the SNP objects
  ## for the rotterdam dataset.
  list_of_objects <- generate_object(exists = exists('list_of_objects'),
                                     dataset = imput_df,
                                     snps = vector_of_snps)

  ## Merging datasets by matching variables names
  matching_variables <- c("ID", "Diag", "Sex", "E4status", "Age_to_use", "Age75", "Centre", vector_of_snps)
  All <- merge(x = df, y = imput_df, by = matching_variables, all.x = T, all.y = T)

  } else {

    cat("Renaming dataset names and subsetting for LOAD patients...\n\n\n")

    ## The only dataset available is the genotyped one
    All <- df
  }

  # Before doing any analysis, making sure the contrasts are performed with respect to the controls, not the Alzheimer cases.
  All$Diag <- relevel(All$Diag, ref = "Control")


  ## Adding binary variables to the merged dataset according to if the data
  ## belong to a specific region, being 0 not belonging and 1 belonging.
  centres <- levels(All$Centre)
  centres <- stringr::str_to_title(centres)
  All[centres] <- 0
  for (centre in centres) {
    All[[centre]][All$Centre == toupper(centre)] <- 1
  }

  if (EP) {
  ## Adding a Region variable that distinguishes between individuals from N.Europe
  ## or from Spain

    if (ADNI) {
      All$Region <- ifelse(All$Centre == "MADRID" | All$Centre == "OVIEDO" | All$Centre == "SANTANDER", "Spain", ifelse(All$Centre == "USA", "America", "N.Eur"))

      } else {
      All$Region <- ifelse(All$Centre == "MADRID" | All$Centre == "OVIEDO" | All$Centre == "SANTANDER", "Spain", "N.Eur")
    }

  } else {
    All$Region <- NA
  }

  ## Obtaining all possible snp inheritance model combinations and adding the new columns to the dataset
  recessive <- paste0(vector_of_snps, "r")
  additive <- paste0(vector_of_snps, "a")
  dominant <- paste0(vector_of_snps, "d")
  inheritance <- c(recessive, additive, dominant)
  All[inheritance] <- NA

  ## Defining the variables (columns) order in the dataset
  col_order <- c("ID", "Diag", "Sex", "Age_to_use", "Age75", "Region",
                 "E4status", "Centre", centres, vector_of_snps, inheritance)

  All <- All[,col_order]

  ## Coding the inheritance with numbers from the genotype values for each snp
  All <- generate_inheritances(snps = list_of_objects, data = All)

  objects_to_be_removed <- c("additive", "dominant", "recessive", "inheritance",
                             "centres", "matching_variables", "col_order",
                             "centre", "snp")
  suppressWarnings(rm(list = objects_to_be_removed))


# 4. Running plink by calling the unix shell -----------------------------

  cat("Running plink from the shell...\n\n\n")

  ## Preliminary steps
  system("mkdir -p Output/plink/0-Data") #Create folder where to store the input
  system('plink --out "$PWD"/Output/plink/0-Data/GPA --file "$PWD"/Input/plink/test --make-bed', ignore.stdout = !verbose) #Make binary files to speed computation
  setwd(file.path(output, "plink"))

  ## Removing missingness
  system('mkdir -p 1-QC/1-Missingness')
  system('plink --bfile 0-Data/GPA --missing --out 1-QC/1-Missingness/missing', ignore.stdout = !verbose)
  system('plink --bfile 0-Data/GPA --geno 0.2 --make-bed --out 1-QC/1-Missingness/GPA_1', ignore.stdout = !verbose)
  system('plink --bfile 1-QC/1-Missingness/GPA_1 --mind 0.2 --make-bed --out 1-QC/1-Missingness/GPA_2', ignore.stdout = !verbose)
  system('plink --bfile 1-QC/1-Missingness/GPA_2 --geno 0.02 --make-bed --out 1-QC/1-Missingness/GPA_3', ignore.stdout = !verbose)
  system('plink --bfile 1-QC/1-Missingness/GPA_3 --mind 0.02 --make-bed --out 1-QC/1-Missingness/GPA_4', ignore.stdout = !verbose)

  ## Removing SNPs with lower MAF than threshold
  system('mkdir 1-QC/3-MAF')
  system('plink --bfile 1-QC/1-Missingness/GPA_4 --maf 0.05 --make-bed --out 1-QC/3-MAF/GPA_7', ignore.stdout = !verbose)


  ## Looking at deviations from Hardy-Weinberg equilibrium
  system('mkdir 1-QC/4-HWE')
  system('plink --bfile 1-QC/3-MAF/GPA_7 --hwe 1e-6 --make-bed --out 1-QC/4-HWE/GPA_8', ignore.stdout = !verbose)
  system('plink --bfile 1-QC/4-HWE/GPA_8 --hwe 1e-10 include-nonctrl --make-bed --out 1-QC/4-HWE/GPA_9', ignore.stdout = !verbose)


  setwd("../..") #Return to the main project directory





# 5. Diagnose possible problems and perform data conversions -------------

  cat("Preparing data for analysis...\n\n\n")

  ## Unexplained level in the E4status variable (coming from Rotterdam dataset)
  All$E4status[All$E4status == 2] <- NA

  ## Convert all categorical variables into the factor class.
  factors <- -c(1,4)
  All[,factors] <- lapply(All[,factors], factor)

  ## Create contingency tables of age, sex and E4status by Dx and Centre
  variables <- c("Age75", "Sex", "E4status", "Diag", "Centre", "Region")

    ### For each combination up to 3 variables at a time get a contingency table
    list_of_tables <- vector(mode = "list", length = 2)

    i = 2 # Two pairs of variables cross tabulation
    while (i < 4) { # While inferior to 4 variables
      # Generate all combinations of variables possible, taking i at a time.
      combinations <- combn(variables, i, simplify = FALSE)
      # Create character vector inside the list of tables to store the names of the tables created
      list_of_tables[[i - 1]] <- character(length = length(combinations))
      # For combination in combinations
      for (n in seq_along(combinations)) {
        # Get the combination
        combination <- combinations[[n]]
        # Create the formula to input to the xtabs function
        formula <- paste0("All$", combination)
        formula <- paste(formula, collapse = " + ")
        formula <- paste0("~", formula)
        # Create the name of the table as variable1_variable2
        name <- paste(combination, collapse = "_")
        # Create the xtab and assign the name created in order to access it.
        assign(name, xtabs(formula))
        # Add the name created to the list of tables to be able to access the table later on.
        list_of_tables[[i - 1]][n] <- name
      }
      # Once all possible 2 variables combinations have been performed, do it
      # for three variables
      i = i + 1
    }
    names(list_of_tables) <- c("Two variables", "Three variables")


  ## Create summary tables of mean, min, max age by Diagnosis and centre
    ### Giving variables subsets names
    ages <- All$Age_to_use
    diagnostic <- All$Diag
    centres <- All$Centre
    regions <- All$Region

    ### Only Age
    variable <- ages
    group <- diagnostic
    AgeAll <- capture.output(psych::describeBy(variable, group, digits = 2))

    ### Age of alzheimer patients by Centre
    variable <- ages[diagnostic == "AD"]
    group <- centres[diagnostic == "AD"]
    AgeAD <- capture.output(psych::describeBy(variable, group, digits = 2))

    ### Age of control patients by Centre
    variable <- ages[diagnostic == "Control"]
    group <- centres[diagnostic == "Control"]
    AgeCTRL <- capture.output(psych::describeBy(variable, group, digits = 2))

    ### Age of alzheimer patients by Region (N.Eur or Spain)
    if (EP) {
    variable <- ages[diagnostic == "AD"]
    group <- regions[diagnostic == "AD"]
    AgeADRegion <- capture.output(psych::describeBy(variable, group, digits = 2))

    variable <- ages[diagnostic == "Control"]
    group <- regions[diagnostic == "Control"]
    AgeCTRLRegion <- capture.output(psych::describeBy(variable, group, digits = 2))
    }

  summary_tables <- c("AgeAll", "AgeAD", "AgeCTRL")
  if (EP) {
    summary_tables <- c(summary_tables, "AgeADRegion", "AgeCTRLRegion")
  }

  for (table_name in summary_tables) {
    text <- c('####################', table_name, '####################', get(table_name), "\n\n")
    cat(text, file = file.path(output, "Age.txt"), sep = "\n", append = TRUE)
  }


  ## Create datasubsets for each region (N.Eur, Spain and UK in the EP dataset)
  if (EP) {
    for (region in levels(All$Region)) {
      subset <- subset(All, Region == region)
      assign(region, subset)
    }
  }


  ## and each centre (Optima, Nottingham, Bristol, Santander, Oviedo, Madrid, Rotterdam)
  for (centre in levels(All$Centre)) {
    subset <- subset(All, Centre == centre)
    centre <- stringr::str_to_title(centre)
    assign(centre, subset)
  }


  ## Now create datasets divided by age, sex and E4status.
  variables <- c("Age75", "Sex", "E4status")
  for (variable in variables) {
    for (level in levels(All[[variable]])) {
      name <- paste(variable, level, sep = "_")
      subset <- All[All[[variable]] == level,]
      assign(name, subset)
    }
  }

  if (EP) {
    ## And the same datasets division but by region
    for (region in levels(All$Region)) {
      Region <- get(region)
      for (variable in variables) {
        for (level in levels(Region[[variable]])) {
          name <- paste(region, variable, level, sep = "_")
          subset <- Region[Region[[variable]] == level,]
          assign(name, subset)
        }
      }
    }
    `Spain_E4status_E4+_Sex_Male` <- subset(Spain, E4status == "E4+" & Sex == "Male")
  }

  ## Subsetting selected centers by variables Age75 or E4status

  if (EP) {
    variables <- variables[-c(2)] #Just keeping Age75 and E4status and dismissing Sex
    centres <- c("Madrid", "Oviedo", "Santander", "Bristol")

    for (centre in centres) {
      Centre <- get(centre)
      for (variable in variables) {
        for (level in levels(All[[variable]])) {
          name <- paste(centre, variable, level, sep = "_")
          subset <- Centre[Centre[[variable]] == level,]
          assign(name, subset)
        }
      }
    }
  }

  ## For each center, subset by variable Sex.

  centres <- str_to_title(levels(All$Centre))
  for (centre in centres) {
    Centre <- get(centre)
    for (level in levels(All$Sex)) {
      name <- paste(centre, level, sep = "_")
      subset <- subset(Centre, Sex == level)
      assign(name, subset)
    }
  }

  for (names in list_of_tables) {
    for (table in names) {
      Table <- get(table)
      if (length(dim(Table)) == 2) {
        crosstab <- as.data.frame.matrix(Table)
      }
      else {
        crosstab <- Table
      }
      write.xlsx2(crosstab, file = file.path(output, "variables_xtabs.xls"), sheetName = table, append = T)
    }
  }

  ## Create xtabs of genotypic counts
  additive_snps <- grep(pattern = "[0-9]a$", colnames(All), value = T) #GET SNPS BY THE ADDITIVE MODEL
  variables <- c("Centre", "Age75")

  if (EP) {
    variables <- c(variables, "Region")
  }

  for (snp in additive_snps) {
    snp <- paste0("All$", snp)
    formula <- paste(snp, "All$Diag", sep = " + ")
    formula <- paste0("~", formula)
    formula <- as.formula(formula)
    if (verbose) {
      cat("\n\n")
      print(sprintf("############# %s #############", snp))
      cat("\n")
      print(sprintf("2-way crosstable of %s by Diagnostic", snp))
      cat("\n")
      print(xtabs(formula))
    }
    for (variable in variables) {
      variable <- paste0("All$", variable)
      formula <- paste("All$Diag", snp, variable, sep = " + ")
      formula <- paste0("~", formula)
      formula <- as.formula(formula)
      value <- xtabs(formula)
      name <- paste(snp, variable, sep = "_")
      assign(name, value)
      if (verbose) {
        cat("\n")
        print(sprintf("3-way crosstable of %s by Diagnostic and %s", snp, variable))
        cat("\n")
        table <- get(name)
        print(ftable(table))
      }
    }
  }

#/////////////////////////////////////ANALYSIS//////////////////////////////////////////////////


# 6. Main Effects Analysis --------
# Establish all main effects, controlling for age, sex, APOE4+/- & study centre
# ,eg Nottingham; apart from APOE4, few main effects, if any, are likely to be
# significant

  cat("Performing main effects analysis...\n\n")

  ## GLM data preparation

    ### Selecting the datasets and variables we want to perform the analysis to
    DATASETS <- "All"
    if (EP) {
      DATASETS <- c(DATASETS, "N.Eur", "Spain")
    }
    if (ADNI) {
      DATASETS <- c(DATASETS, "America")
    }
    variables <- c("Sex", "E4status", "Age75", "Centre")

    ### Order datasets variables by alfabetical name of the variables and replace the original ones
    for (i in seq_along(DATASETS)) {
      assign(DATASETS[i], get(DATASETS[i])[,order(names(get(DATASETS[i])))])
    }

    if (EP) {
      ### Relevel Santander instead of Madrid as the reference group for the glm as Santander has lower population n
      Spain$Centre <- relevel(Spain$Centre, ref = 6)
    }

    master_list <- perform_analysis(.mode = "main_effects", .data = DATASETS, snps = list_of_objects, covariates = variables)


  ## Print and count significant results obtained for each dataset
   if (verbose) {
      print_significant_results(master_list)
    }


# 7. Variable-SNP interaction analysis -------------------------------------------------


# Examine all possible interactions with age 75y, sex & APOE4;
# this means repeating the models from (i) adding one interaction term at a time.
# Since synergy factors work with binary, not continuous, variables, I've specified age 75y, sex & APOE4

  cat("\n\nPerforming snp-variable interaction analysis...\n\n")

  master_list <- perform_analysis(.mode = "interaction", .data = DATASETS, snps = list_of_objects, covariates = variables)

  list_of_covariate_snp_pairs <- vector(mode = "list")

  ## Extract significant results in the All dataset
  for (i in seq_along(master_list)) {
    dataset_list <- master_list[[i]]
    if (names(master_list)[i] == "All") {
      for (j in seq_along(dataset_list)) {
        snp_list <- dataset_list[[j]]
        for (n in seq_along(snp_list[["Interactions"]][["Other_covariates"]])) {
          covariate_list <- snp_list[["Interactions"]][["Other_covariates"]][[n]]
          for (x in seq_along(covariate_list[["SF"]])) {
            if (covariate_list[["Significant"]][x] == T) {
              if (names(master_list)[i] == "All") {
                ma_list <- list(names(master_list)[i],
                                names(dataset_list)[j],
                                snp_list[["Gene"]],
                                names(snp_list[["Interactions"]][["Other_covariates"]])[n],
                                covariate_list[["Name"]],
                                covariate_list[["SF"]],
                                covariate_list[["Significant"]],
                                tail(covariate_list[["Summary"]]),
                                tail(snp_list[["Main_effects"]]))
                if (verbose) {
                  lapply(ma_list, print)
                  print("--------------------------------")
                }
            list_of_covariate_snp_pairs <- append(x = list_of_covariate_snp_pairs, values = paste(names(dataset_list)[j], names(snp_list[["Interactions"]][["Other_covariates"]])[n], sep = "_"))
            }
          }
        }
      }
    }
  }
}


# 7.5 Examine subsets of found interactions.
  if (verbose) {
    cat("Examining subsets of found interactions...\n\n")
  for (item in list_of_covariate_snp_pairs) {
    snp <- unlist(strsplit(item, split = "_"))[1]
    covariate <- unlist(strsplit(item, split = "_"))[2]
    covariates <- setdiff(variables, covariate)
    levels <- levels(All[[covariate]])
    datasets <- paste(covariate, levels, sep = "_")
    for (dataset in datasets) {
      perform_subset_analysis(dataset = dataset, snp = snp, covariate = covariate, covariates = covariates, verbose = verbose)
      if (verbose) {
        cat("\n#### Interaction results ####\n")
        print(tail(master_list$All[[snp]]$Interactions$Other_covariates[[covariate]]$Summary, n = 5))
      }
    }
  }
}


# 8. SNP-SNP interaction analysis --------------------------------------------

# Examine all possible interactions between snps from a selected list of snp-snp pairs
#
  cat("\n\nPerforming snp-snp interaction analysis...\n\n")

  snp_snp <- read.csv(file = file.path(input, "snp-snp_selection.csv"))
  snp_snp <- snp_snp[-c(30:32),]
  rownames(snp_snp) <- NULL

  ## Analyze SNP-SNP interactions with reported models
  list_of_snp_pairs <- generate_list_of_snp_pairs(mode = "reported_model",
                                                  snp_df = snp_snp)

  master_list <- store_interactions(snp_dataset = snp_snp, master_list = master_list)

  # Performing the analysis
  master_list <- perform_analysis(.mode = "snp",
                                  .submode = "reported_model",
                                  .data = DATASETS,
                                  snps = list_of_snp_pairs,
                                  covariates = variables,
                                  .verbose = T)

  # PRINT STEP

  significant_snp_pairs <- vector(mode = "character")

  if (verbose) {
  for (i in seq_along(master_list)) {
    cat("\n\n\n", "#########################", names(master_list)[i], "#########################", "\n\n\n")
    for (snp_pair in list_of_snp_pairs) {
      snp1_model <- snp_pair[1]
      snp1 <- sub(snp1_model, pattern = "[a-z]$", replacement = "")
      snp2_model <- snp_pair[2]
      snp2 <- sub(snp2_model, pattern = "[a-z]$", replacement = "")
      significance <- master_list[[i]][[snp1]][["Interactions"]][["SNP"]][["Reported"]][[snp1_model]][[snp2]][[snp2_model]][["Significant"]]
      for (x in seq_along(significance)) {
        if (significance[x] == TRUE) {
          cat("---", snp1, "---", "\n\n")
          gene <- master_list[[i]][[snp1]][["Gene"]]
          covariate <- master_list[[i]][[snp2]][["Gene"]]
          print(sprintf("Genes are %s and %s", gene, covariate))
          result <- master_list[[i]][[snp1]][["Interactions"]][["SNP"]][["Reported"]][[snp1_model]][[snp2]][[snp2_model]]
          print(result)
          if (names(master_list)[i] == "All") {
            significant_snp_pairs <- unique(append(significant_snp_pairs, result$Name[x]))
          }
        }
      }
    }
  }
}

  ## Analyze SNP-SNP interactions with best models
  list_of_snp_pairs2 <- generate_list_of_snp_pairs(mode = "best_model",
                                                  snp_df = snp_snp)

    ## Creating the containers that will store the interactions from the analysis
    master_list <- store_possible_interactions(.master_list = master_list, pairs = list_of_snp_pairs2)

    ## Performing the analysis
    master_list <- perform_analysis(.mode = "snp",
                                    .submode = "best_model",
                                    .data = DATASETS,
                                    snps = list_of_snp_pairs2,
                                    covariates = variables,)

    # PRINT STEP

    significant_snp_pairs2 <- vector(mode = "character")

    if (verbose) {
      for (i in seq_along(master_list)) {
        cat("\n\n\n", "#########################", names(master_list)[i], "#########################", "\n\n\n")
        for (snp_pair in list_of_snp_pairs2) {
          snp1 <- snp_pair[1]
          snp2 <- snp_pair[2]
          significance <- master_list[[i]][[snp1]][["Interactions"]][["SNP"]][["Best_model"]][[snp2]][["Significant"]]
          for (x in seq_along(significance)) {
            if (significance[x] == TRUE) {
              cat("---", snp1, "---", "\n\n")
              gene <- master_list[[i]][[snp1]][["Gene"]]
              covariate <- master_list[[i]][[snp2]][["Gene"]]
              print(sprintf("Genes are %s and %s", gene, covariate))
              result <- master_list[[i]][[snp1]][["Interactions"]][["SNP"]][["Best_model"]][[snp2]]
              print(result)
              if (names(master_list)[i] == "All") {
                significant_snp_pairs2 <- unique(append(significant_snp_pairs2, result$Name[x]))
              }
            }
          }
        }
      }
    }

# 9. Plotting -------------------------------------------------------------

  if (verbose) {

  # 9.1 MAIN EFFECTS

  plot_df <- data.frame(matrix(NA, nrow = length(list_of_objects), ncol = 5))
  colnames(plot_df) <- c("Gene", "pvalue", "OR", "diff_expressed", "labeled")

  for (i in seq_along(master_list$All)) {
    gene <- master_list$All[[i]]$Gene
    pvalue <- unname(unlist(tail(master_list$All[[i]]$Main_effects, n = 1)[3]))
    OR <- unname(unlist(tail(master_list$All[[i]]$Main_effects, n = 1)[2]))
    diff_expressed <- ifelse(pvalue < 0.05 & {OR > 1.2 | OR < 0.8}, TRUE, FALSE)
    labeled <- ifelse(diff_expressed == TRUE, gene, NA)
    plot_df[i,1] <- gene
    plot_df[i,2] <- pvalue
    plot_df[i,3] <- OR
    plot_df[i,4] <- diff_expressed
    plot_df[i,5] <- labeled
  }
  p <- ggplot2::ggplot(data = plot_df, aes(x = OR, y = -log10(pvalue), col = diff_expressed, label = labeled)) + geom_point() + theme_minimal() + geom_text(vjust = -1) + coord_cartesian(xlim = c(0, 9), ylim = c(0, 1000)) + guides(color = FALSE)

  #p2 <- p + geom_vline(xintercept = c(0.8, 1.2), col = "red") + geom_hline(yintercept = -log10(0.05), col = "red")
  p
  }
  # 9.2 SNP-VARIABLE INTERACTION

    ## 9.2.1 Extracting results into a dataframe
    plot_df <- data.frame(matrix(NA, nrow = length(list_of_covariate_snp_pairs), ncol = 7))
    colnames(plot_df) <- c("snp", "SF", "lower", "upper", "p_value", "gene", "covariate")
    for (pair in list_of_covariate_snp_pairs) {
      i <- which(list_of_covariate_snp_pairs == pair)
      snp <- strsplit2vector(pair, "_")[1]
      covariate <- strsplit2vector(pair, "_")[2]
      summary <- master_list$All[[snp]]$Interactions$Other_covariates[[covariate]]$Summary
      if (grepl("a", names(master_list$All[[snp]]$Best_model))) {
        summary <- tail(summary, n = 2)
      } else {
        summary <- tail(summary, n = 1)
      }
      result <- summary[summary[,6] < 0.05]
      plot_df[i, 1] <- snp
      plot_df[i, 2] <- result[3]
      plot_df[i, 3] <- result[4]
      plot_df[i, 4] <- result[5]
      plot_df[i, 5] <- result[6]
      plot_df[i, 6] <- master_list$All[[snp]]$Gene
      plot_df[i, 7] <- covariate
    }
    p <- ggplot(data = plot_df, aes(x = SF, y = snp, xmin = lower, xmax = upper))
    p <- p + geom_pointrange(aes(col = gene))
    p <- p + geom_vline(aes(fill = gene), xintercept = 1, linetype = "dotted")
    p <- p + geom_errorbar(aes(xmin = lower, xmax = upper, col = gene),width = 0.5, cex = 1)
    p <- p + facet_wrap(~covariate, strip.position = "left", nrow = 3, scales = "free_y")
    p <- p + geom_text(aes(label = round(p_value, 3)),hjust = -0.5, vjust = -1)
    p <- p + scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4))

  # 9.3 SNP-SNP INTERACTION

    ## 9.3.1 Extracting results into a dataframe
    plot_df <- data.frame(matrix(NA, nrow = length(significant_snp_pairs)*length(master_list), ncol = 11))
    colnames(plot_df) <- c("snp1", "gene1", "snp2", "gene2", "SF", "lower", "upper", "p_value", "interaction", "dataset", "gene_interaction")
    for (interaction in significant_snp_pairs) {
      position <- which(significant_snp_pairs == interaction)
      snp_pair <- strsplit2vector(interaction, ":")
      snp_1 <- sub(pattern = "[a-z][0-9]$", "", snp_pair[1])
      snp_2 <- sub(pattern = "[a-z][0-9]$", "", snp_pair[2])
      gene1 <- master_list$All[[snp_1]]$Gene
      gene2 <- master_list$All[[snp_2]]$Gene
      for (i in seq_along(master_list)) {
        interaction <- master_list[[i]][[snp_1]]$Interactions$SNP[[snp_2]]
        summary <- interaction$Summary
        if (grepl("a", interaction$Name)) {
          summary <- tail(summary, n = 2)
        } else {
          summary <- tail(summary, n = 1)
        }
        index <- (i - 1) * length(significant_snp_pairs) + position
        if (i == 1) {
          result <- summary[summary[,6] < 0.05]
        } else {
          result <- summary[1,]
        }
        plot_df[index, 1] <- snp_1
        plot_df[index, 2] <- gene1
        plot_df[index, 3] <- snp_2
        plot_df[index, 4] <- gene2
        plot_df[index, 5] <- result[3]
        plot_df[index, 6] <- result[4]
        plot_df[index, 7] <- result[5]
        plot_df[index, 8] <- result[6]
        plot_df[index, 9] <- paste(snp_1, snp_2, sep = "*")
        plot_df[index, 10] <- names(master_list)[i]
        plot_df[index, 11] <- paste(gene1, gene2, sep = "*")
      }
    }

    p <- ggplot(data = plot_df, aes(x = dataset, y = SF, ymin = lower, ymax = upper))
    p <- p + geom_pointrange(aes(col = dataset))
    p <- p + geom_hline(aes(fill = dataset), yintercept = 1, linetype = "dotted")
    p <- p + geom_errorbar(aes(ymin = lower, ymax = upper, col = dataset), width = 0.5, cex = 1)
    p <- p + facet_wrap(~gene_interaction, strip.position = "left", nrow = length(significant_snp_pairs), scales = "free_y")
    p <- p + geom_text(aes(label = round(p_value, 3)),hjust = 0, vjust = -1)
    p <- p + coord_flip()
    p <- p+theme(axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
    p <- p + scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4))

    vector_of_sf_replicas <- vector(mode = "character", length = nrow(snp_snp))
    vector_of_p_values <- vector(mode = "numeric", length = nrow(snp_snp))
    # TABLE OF RESULTS FOR SELECTED INTERACTIONS
    for (i in seq_len(nrow(snp_snp))) {
      row <- snp_snp[i,]
      snp1 <- row$Genotype_1
      snp2 <- row$Genotype_2
      summary <- master_list$All[[snp1]]$Interaction$SNP[[snp2]]$Summary
      result <- tail(summary, n = 1)
      vector_of_sf_replicas[i] <- round(result[3], digits = 2)
      vector_of_p_values[i] <- round(result[6], digits = 2)
    }

  snp_snp_results <- snp_snp
  snp_snp_results$SF_replica <- vector_of_sf_replicas
  snp_snp_results$p_values_replica <- vector_of_p_values

  snp_snp_results$Comments <- sub(" \\(.*$", "", x = snp_snp_results$Comments)
  cat("\n\nDONE!\n\n")

  vector_of_sf_replicas <- vector(mode = "character", length = nrow(snp_apoe))
  vector_of_p_values <- vector(mode = "numeric", length = nrow(snp_apoe))
  for (i in seq_len(nrow(snp_apoe))) {
    row <- snp_apoe[i,]
    snp1 <- row$Genotype_1
    summary <- master_list$All[[snp1]]$Interaction$Other_covariates$E4status$Summary
    result <- tail(summary, n = 1)
    vector_of_sf_replicas[i] <- round(result[3], digits = 2)
    vector_of_p_values[i] <- round(result[6], digits = 2)
  }
  snp_apoe_results <- snp_apoe
  snp_apoe_results$SF_replica <- vector_of_sf_replicas
  snp_apoe_results$p_values_replica <- vector_of_p_values

  CDK5_forest_plot <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 6))
  colnames(CDK5_forest_plot) <- c("SF", "lower", "upper", "p_value", "dataset", "gene")
  for (i in seq_along(master_list)) {
    results <- master_list[[i]]$rs2069442$Interactions$Other_covariates$E4status$Summary
    results <- tail(results, n=1)
    CDK5_forest_plot[i, 1] <- results[3]
    CDK5_forest_plot[i, 2] <- results[4]
    CDK5_forest_plot[i, 3] <- results[5]
    CDK5_forest_plot[i, 4] <- results[6]
    CDK5_forest_plot[i, 5] <- names(master_list)[i]
    CDK5_forest_plot[i, 6] <- "CDK5"
  }

  p <- ggplot(data = CDK5_forest_plot, aes(x = SF, y = dataset, xmin = lower, xmax = upper))
  p <- p + geom_pointrange(aes(col = dataset))
  p <- p + geom_vline(aes(fill = dataset), xintercept = 1, linetype = "dotted")
  p <- p + geom_errorbar(aes(xmin = lower, xmax = upper, col = dataset),width = 0.5, cex = 1)
  p <- p + facet_wrap(~gene, strip.position = "left", nrow = 1, scales = "free_y")
  p <- p + geom_text(aes(label = round(p_value, 3)),hjust = -0.5, vjust = -1)
  p <- p+theme(axis.title.y=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
  p <- p + scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4))

