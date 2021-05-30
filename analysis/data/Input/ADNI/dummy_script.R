
data <- c("duplicated_a1", "duplicated_a2", "duplicated_a3", "duplicated_a4")

list_of_duplications <- vector(mode = "list", length = length(duplicated_iid))
names(list_of_duplications) <- duplicated_iid
for (dupli in duplicated_iid) {
  snp_list <- vector(mode = "list", length = length(snp_id))
  names(snp_list) <- snp_id
  for (snp in snp_id) {
    dataset_list <- vector("list", length(data))
    names(dataset_list) <- data
    for (df in data) {
      name <- df
      df <- get(df)
      if (dupli %in% df$iid) {
        if (snp %in% colnames(df[7:ncol(df)])) {
          value <- df[[snp]][df$iid == dupli]
          dataset_list[[name]] <- value
        }
      }
    }
    snp_list[[snp]] <- dataset_list
  }
  list_of_duplications[[dupli]] <- snp_list
}

for (dupli in names(list_of_duplications)) {
  for (snp in names(list_of_duplications[[dupli]])) {
    inner_list <- list_of_duplications[[dupli]][[snp]]
    list_of_duplications[[dupli]][[snp]] <- unlist(inner_list)
  }
}

for (dupli in names(list_of_duplications)) {
  for (snp in names(list_of_duplications[[dupli]])) {
    table <- table(list_of_duplications[[dupli]][[snp]])
    if (length(table) != 1) {
      print(dupli)
      print(snp)
      print(table)
      print("------")
    }
  }
}
