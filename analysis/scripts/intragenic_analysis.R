## 8.5 In between snp analysis

  vector_of_genes <- sapply(list_of_objects, function(x) x$gene)
  vector_of_genes <- unique(unname(vector_of_genes))

  list_of_inbetween_snps <- vector(mode = "list", length = length(vector_of_genes))
  for (i in seq_along(vector_of_genes)) {
    gene <- vector_of_genes[i]
    snps <- unlist(lapply(list_of_objects, function(x) if (x$gene == gene) {x$id}))
    if (length(snps) > 1) {
      list_of_inbetween_snps[[i]] <- snps
      names(list_of_inbetween_snps)[i] <- gene
    } else {
      next
    }
  }
  not_null <- unname(sapply(list_of_inbetween_snps, function(x) !is.null(x)))
  list_of_inbetween_snps <- list_of_inbetween_snps[not_null]
  vector_of_combinations <- vector(mode = "character", length = length(list_of_inbetween_snps))

  for (i in seq_along(list_of_inbetween_snps)) {
    gene <- names(list_of_inbetween_snps)[i]
    combinations <- combn(list_of_inbetween_snps[[i]], m = 2, simplify = F)
    combinations <- lapply(combinations, unname)
    name <- paste(gene, "combinations", sep = "_")
    assign(name, combinations)
    vector_of_combinations[i] <- name
  }

  for (snp_combinations in vector_of_combinations) {
    list_snp_combinations <- get(snp_combinations)
    list_of_snp_pairs <- c(list_of_snp_pairs, list_snp_combinations)
  }