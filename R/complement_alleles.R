#' Get the complementary allele sequence
#'
#' For a character vector of a genomic sequences it returns a character vector
#' with the complementary sequence.
#'
#' @param sequence character vector, with a nucleotide sequence.
#' @param mode string. If "dna" it complements "A" to "T". If "rna" it
#' complements "A" to "U".
#'
#' @return character vector
#' @export
complement_alleles <-
  function(sequence, mode = "dna") {

    # Creating complement function
    complement_function <- function(allele, .mode = mode) {
      if (allele == "G") {
        complement <- "C"
      }
      if (allele == "C") {
        complement <- "G"
      }
      if (allele == "A") {
        if (mode == "dna") {
          complement <- "T"
        } else if (mode == "rna") {
          complement <- "U"
        }
      }
      if (allele == "T" | allele == "U") {
        complement <- "A"
      }
      complement
    }
    # Using complement function on sequence elements
    sequence <- toupper(sequence)
    complementary <- sapply(sequence, complement_function)
    complementary <- unname(complementary)

    # Returning complementary sequence obtained
    complementary
  }
