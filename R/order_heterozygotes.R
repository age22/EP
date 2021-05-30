#' Order heterozygotes by alphabetical order
#'
#' @param x character vector, consisting on the genotypes we want to order
#' alphabetically
#'
#' @return character vector, ordered list of genotypes
#' (the heterozygotes are ordered alphabetically)
#' @export
order_heterozygotes <- function(x) {
  sort_genotypes <- function(y) {
    if (is.na(y)) {
      genotype <- NA
    } else {
      alleles <- strsplit2vector(y, pattern = "")
      sorted_alleles <- sort(alleles)
      genotype <-  paste(sorted_alleles, collapse = "")
    }
    genotype
  }
  sapply(x, sort_genotypes)
}



