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
