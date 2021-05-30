#' Recodifies an imputated dataframe into genotypes
#'
#' For a given imputated dataframe, it recodifies the output obtained into
#' genotype data according to a reference panel and checks that it matches the
#' same strand in the genotyped dataframe through the list_of_objects.
#'
#' @param list_of_objects 'SNP_set' S3 class, contains several SNP objects with
#' their respective major and minor allele that will be used to see if they
#' match the reference and alternative allele in the reference panel
#' (thousand genomes)
#' @param df dataframe, with the imputated data for each snp under study
#' @param thousand_genomes data.frame, which contains the reference panel
#' reference and alternative allele with respect to which the imputation took
#' place
#' @param verbose logical, should the function print to output the steps it
#' performs?
#' @param ... further arguments passed to or from other methods.
#' @return A dataframe with the imputation converted into genotypes.
#' @export
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
