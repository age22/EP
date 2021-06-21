#' Creates a list of SNP objects
#'
#' Given a data frame, it checks if a list of SNP objects exist and if not
#' creates it with a set of given snp id's. If it does exist it expects that
#' the data frame inputted is an imputated one and for each SNP adds its
#' attributes to the SNP object.
#'
#' @param exists logical, does the list of SNP objects already exist?
#' @param snps character vector with "gene_snp" names in the case of genotyped
#' df (used to retrieve the gene names) or just RefSNP id's if imputated.
#' @param ... further arguments passed to or from other methods.
#' Use argument dataset, for inputting a dataframe, with all the SNP data and
#' argument verbose, logical, for whether the function should print output.
#'
#' @return a SNP_set S3 object (a list of SNP objects)
#' @export
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
  function(snp, dataset, ...) {
    ## Generating genotypic counts by doing a summary table of the
    ## genotypes in the control cases for a specific snp
    controls_subset <- dataset[dataset$Diag == "Control",]
    snp_controls <- controls_subset[[snp]]
    geno_count <- table(snp_controls)

    ## Calculate genotype frequencies
    geno_freq <- 100*geno_count/sum(geno_count)
    geno_freq <- round(geno_freq, digits = 3)
    ## Calculate allele counts
    allele_count <- calculate_allele_counts(geno_count, dataset, snp, ...)

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
  function(genotype_count, dataset, snp, verbose) {
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


#' @export
print.SNP <-
  function(x, ...) {
    cat("######", x$id, "######")
    for (i in seq_along(x)) {
      attribute <- names(x)[i]
      value <- x[[i]]
      name <- names(x[[i]])
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
    invisible(NULL)
  }

#' @export
print.SNP_set <-
  function(x, ...) {
    cat("This is just a convenience function, invocated when print() is used, to provide an overview of the genes and snps present in this list. To view its contents, use list_of_objects$your_snp_id\n\n\n")
    snps <- names(x)
    genes <- unique(unname(sapply(x, function(y) y$gene)))
    print(list('snps' = snps, 'genes' = genes))
    cat(sprintf("\nThe list has a total of %s snps and %s genes", length(snps), length(genes)))
    invisible(NULL)
  }

#' @export
str.SNP_set <-
  function(object, ...) {
    cat("This is just a convenience function, invocated when str() is used, to provide an overview of the list elements. To view its contents, use list_of_objects$your_snp_id\n\n")
    cat("An example of the structure of the objects in the list is the following one:\n\n")
    print(utils::ls.str(object[1]))
    invisible(NULL)
  }
