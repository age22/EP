# MAIN EFFECTS

snp_pvalue_list <- vector(mode = "list", length = length(master_list$All))
for (snp in names(master_list$All)) {
  Gene <- list_of_objects[[snp]]$gene
  best_model <- names(master_list$All[[snp]]$Best_model)
  p_values <- master_list$All[[snp]]$Main_effects$Pr...t..
  p_value <- tail(p_values, n = 1)
  i <- which(names(master_list$All) == snp)
  snp_pvalue_list[[i]] <- p_value
  names(snp_pvalue_list)[i] <- snp
}

p_values <- unlist(snp_pvalue_list, use.names = F)
p_adjusted <- p.adjust(p_values, method = "BH")
for (snp in names(snp_pvalue_list[c(which(p_adjusted < 0.05))])) {
  print(snp)
  print(list_of_objects[[snp]]$gene)
  cat("\n")
}

# SNP-SNP INTERACTION

## Best model
library(rlist)
snp_pvalue_list <- vector(mode = "list", length = length(list_of_snp_pairs2))
for (i in seq_along(list_of_snp_pairs2)) {
  pair <- list_of_snp_pairs2[[i]]
  snp1 <- pair[1]
  snp2 <- pair[2]
  interaction <- master_list$All[[snp1]]$Interactions$SNP$Best_model[[snp2]]
  interaction_name <- interaction$Name
  if (is.null(interaction$Summary[[1]])) {
    next
  }
  p_values <- unname(interaction$Summary[,6])
  p_value <- tail(p_values, n = 1)
  snp_pvalue_list[[i]] <- p_value
  names(snp_pvalue_list)[i] <- interaction_name[1]
}
snp_pvalue_list2 <- rlist::list.clean(snp_pvalue_list)
p_values2 <- unlist(snp_pvalue_list, use.names = F)
p
## Reported models
snp_pvalue_list <- vector(mode = "list", length = length(list_of_snp_pairs))
for (i in seq_along(list_of_snp_pairs)) {
  pair <- list_of_snp_pairs[[i]]
  snp1_model <- pair[1]
  snp2_model <- pair[2]
  snp1 <- sub(snp1_model, pattern = "[a-z]$", replacement = "")
  snp2 <- sub(snp2_model, pattern = "[a-z]$", replacement = "")
  interaction <- master_list$All[[snp1]]$Interactions$SNP$Reported[[snp1_model]][[snp2]][[snp2_model]]
  interaction_name <- interaction$Name
  if (is.null(interaction$Summary[[1]])) {
    next
  }
  p_values <- unname(interaction$Summary[,6])
  p_value <- tail(p_values, n = 1)
  snp_pvalue_list[[i]] <- p_value
  names(snp_pvalue_list)[i] <- interaction_name[1]
}
p_values <- unlist(snp_pvalue_list, use.names = F)
p_adjusted <- p.adjust(p_values, method = "BH")

## SNP-COVARIATE INTERACTION

snp_pvalue_list <- vector(mode = "list", length = length(master_list$All))
covariate <- c("E4status")
for (snp in names(master_list$All)) {
    interaction <- master_list$All[[snp]]$Interactions$Other_covariates[[covariate]]
    p_values <- interaction$Summary[,6]
    p_value <- tail(p_values, n = 1)
    i <- which(names(master_list$All) == snp)
    snp_pvalue_list[[i]] <- p_value
    names(snp_pvalue_list)[i] <- paste(snp, covariate, sep = "_")
  }


p_values <- unlist(snp_pvalue_list, use.names = F)
p_adjusted <- p.adjust(p_values, method = "BH")
for (snp in names(snp_pvalue_list[c(which(p_adjusted < 0.05))])) {
  print(snp)
  print(list_of_objects[[snp]]$gene)
  cat("\n")
}