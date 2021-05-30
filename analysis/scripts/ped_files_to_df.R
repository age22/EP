alleles_start <- 7
number_of_cols <- ncol(my_df)
new_col <- number_of_cols + 1
for (i in alleles_start:number_of_cols) {
  if (i %% 2 == 1) {
    genotypes <- paste0(my_df[, i], my_df[, i+1])
    my_df[new_col] <- genotypes
    new_col <- new_col + 1
  }
}

my_df <- my_df[,-c(alleles_start:number_of_cols)]