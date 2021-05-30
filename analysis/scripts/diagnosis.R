## CHECK MA CONCORDANCE BETWEEN GENOTYPED AND IMPUTED AND IF NOT OUTPUT ALLELE FREQUENCIES ##
wrong <- vector(mode = "list", length = length(list_of_objects))
counter = 0
for (i in seq_along(list_of_objects)) {
  cat("\n\n--------------------", list_of_objects[[i]]$id, "--------------------\n")
  if (list_of_objects[[i]]$minor_allele == list_of_objects[[i]]$RT_minor_allele) {
    cat("AGREEMENT\n")
    if (list_of_objects[[i]]$id %in% weird) {
      print(list_of_objects[[i]]$minor_allele)
      print(list_of_objects[[i]]$RT_minor_allele)
      counter = counter + 1
      cat("\nGenotyped:\n\n")
      print(list_of_objects[[i]]$allele_freq)
      cat("\n\nImputated:\n\n")
      print(list_of_objects[[i]]$RT_allele_freq)
      }
  } else {
    cat("DISAGREEMENT")
    cat("\nGenotyped:\n\n")
    print(list_of_objects[[i]]$allele_freq)
    cat("\n\nImputated:\n\n")
    print(list_of_objects[[i]]$RT_allele_freq)
    wrong[[i]] <- list_of_objects[[i]]$id
  }
}
print(sprintf("There are %s snps in disagreeement in allele frequencies", counter))
wrong <- unlist(wrong)
wrong <- wrong[is.null(wrong) == FALSE]

raro <- c("rs670139", "rs709149")

## DIFFERENCE IN SNPS IN GENOTYPED DATAFRAME AND THE IMPUTATED ONE
setdiff(vector_of_snps, colnames(df2))

# Results = the snp rs70991108, belonging to the gene DHFR doesn't exist in the 
# Rotterdam dataset