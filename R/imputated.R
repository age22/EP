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
#' @param scheme data.frame, which contains the reference panel
#' reference and alternative allele with respect to which the imputation took
#' place
#' @param match_strands boolean, should the program check to match strands with
#' list_of_objects existing alleles?
#' @param verbose boolean, should the function print to output the steps it
#' performs?
#' @return A dataframe with the imputation converted into genotypes.
#' @export
genotype_imputated_df <-
  function(list_of_objects, df, scheme, match_strands = F, verbose = F) {

    if (match_strands) {
      # Logical vector giving info about for each SNP whether it is possible to
      # match two sets of alleles, having each set a pair of alleles. If its
      # possible to match the sets directly or through the complementary alleles
      # the value will be TRUE, else it will be FALSE.
      #
      # Example1 set1 = ("A", "G"), set2 = ("A", "G") is TRUE
      # Example2 set1 = ("A", "G"), set2 = ("C", "T") is TRUE
      # Example3 set1 = ("A", "G"), set2 = ("T", "G") is FALSE
      match_possible <- vector(mode = "list", length = length(list_of_objects))
    }

    for (i in seq_along(list_of_objects)) {
      SNP <- list_of_objects[[i]]
      snp <- SNP$id

      index <- scheme$snp_id == snp
      ref_allele <- scheme[index,]$reference
      alt_allele <- scheme[index,]$alternative
      major_allele <- SNP$major_allele
      minor_allele <- SNP$minor_allele
      genotyped_alleles <- c(major_allele, minor_allele)

      if (match_strands) {
        if (verbose) {
          print(sprintf("Confirming that there is agreement in the chromosomal strand between the imputed data and the genotyped data for the snp %s...", snp))
        }

        scheme_alleles <- c(ref_allele, alt_allele)

        set1 <- genotyped_alleles
        set2 <- scheme_alleles

        match_possible[[i]] <- check_strand_concordance(set1, set2)
        names(match_possible)[i] <- snp
      }

      ### Allele dosage data is recoded into one of the three genotypes possible
      ### in which a value higher than 1.5 means that the individual is
      ### homozygote for allele 1, if lower than 0.5 homozygote for allele 2
      ### and if in between heterozygote


      allele_1 <- ifelse(ref_allele %in% genotyped_alleles,
                         ref_allele,
                         complement_alleles(ref_allele))

      allele_2 <- ifelse(alt_allele %in% genotyped_alleles,
                         alt_allele,
                         complement_alleles(alt_allele))

      homozygote_1 <- paste0(allele_1, allele_1)
      homozygote_2 <- paste0(allele_2, allele_2)
      heterozygote <- paste0(allele_1, allele_2)

      df[[snp]][df[[snp]] > 1.5] <- homozygote_1
      df[[snp]][df[[snp]] < 0.5] <- homozygote_2
      df[[snp]][df[[snp]] >= 0.5 & df[[snp]] <= 1.5] <- heterozygote
    }
    df
  }

#' Check strand concordance
#'
#'
#' @param set1 with two alleles
#' @param set2 with another two alleles
#'
#' @return boolean, TRUE or FALSE depending on whether
#' @export
check_strand_concordance <-
  function(set1, set2) {
    set3 <- complement_alleles(set2)
    if (setequal(set1,set2) == FALSE & setequal(set1, set3) == FALSE) {
      match_possible <- FALSE
    } else {
      match_possible <- TRUE
    }
  }
