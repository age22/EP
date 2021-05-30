CYP46A1_interaction <- vector(mode = "list", length = length(levels(All$Centre)))
for (centre in str_to_title(levels(All$Centre))) {
  snp1_model <- "rs7157609r"
  snp2_model <- "rs4900442r"
  covariates <- c("Sex", "Age75", "E4status")
  interaction <- paste0(snp1_model, "*", snp2_model)
  predictors <- c(covariates, interaction)
  formula <- as.formula(paste("Diag", paste(predictors, collapse = " + "), sep = " ~ "))
  args <-  list(formula = formula, family = quasibinomial("logit"), data = as.name(centre))
  linear_model <- do.call(glm, args)
  OR_summary <- Coeff.OR2(linear_model)
  cat("\n")
  print(centre)
  print(OR_summary)
  i <- which(str_to_title(levels(All$Centre)) == centre)
  CYP46A1_interaction[[i]] <- OR_summary
  names(CYP46A1_interaction)[i] <- centre
}
 #PLOT 

forest_plot <- as.data.frame(matrix(data = NA, nrow = length(CYP46A1_interaction), ncol = 6))
colnames(forest_plot) <- c("SF", "lower", "upper", "p_value", "Centre", "gene")
for (i in seq_along(CYP46A1_interaction)) {
  results <- CYP46A1_interaction[[i]]
  results <- tail(results, n=1)
  forest_plot[i, 1] <- results[3]
  forest_plot[i, 2] <- results[4]
  forest_plot[i, 3] <- results[5]
  forest_plot[i, 4] <- round(results[6], digits = 2)
  forest_plot[i, 5] <- names(CYP46A1_interaction)[i]
  forest_plot[i, 6] <- "CYP46A1*CYP46A1"
}
p <- ggplot(data = forest_plot, aes(x = SF, y = Centre, xmin = lower, xmax = upper))
p <- p + geom_pointrange(aes(col = Centre))
p <- p + geom_vline(aes(fill = Centre), xintercept = 1, linetype = "dotted")
p <- p + geom_errorbar(aes(xmin = lower, xmax = upper, col = Centre),width = 0.5, cex = 1)
p <- p + facet_wrap(~gene, strip.position = "left", nrow = 1, scales = "free_y")
p <- p + geom_text(aes(label = round(p_value, 3)),hjust = -0.5, vjust = -1)
p <- p+theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank())
p



APOE4_CDK5_interaction <- vector(mode = "list", length = length(levels(All$Centre)))
for (centre in str_to_title(levels(All$Centre))) {
  snp1_model <- "rs2069442r"
  snp2_model <- "E4status"
  covariates <- c("Sex", "Age75")
  interaction <- paste0(snp1_model, "*", snp2_model)
  predictors <- c(covariates, interaction)
  formula <- as.formula(paste("Diag", paste(predictors, collapse = " + "), sep = " ~ "))
  args <-  list(formula = formula, family = quasibinomial("logit"), data = as.name(centre))
  linear_model <- do.call(glm, args)
  OR_summary <- Coeff.OR2(linear_model)
  cat("\n")
  print(centre)
  print(OR_summary)
  i <- which(str_to_title(levels(All$Centre)) == centre)
  APOE4_CDK5_interaction[[i]] <- OR_summary
  names(APOE4_CDK5_interaction)[i] <- centre
}

#plot
forest_plot <- as.data.frame(matrix(data = NA, nrow = length(APOE4_CDK5_interaction), ncol = 6))
colnames(forest_plot) <- c("SF", "lower", "upper", "p_value", "Centre", "gene")
for (i in seq_along(APOE4_CDK5_interaction)) {
  results <- APOE4_CDK5_interaction[[i]]
  results <- tail(results, n=1)
  forest_plot[i, 1] <- results[3]
  forest_plot[i, 2] <- results[4]
  forest_plot[i, 3] <- results[5]
  forest_plot[i, 4] <- round(results[6], digits = 2)
  forest_plot[i, 5] <- names(APOE4_CDK5_interaction)[i]
  forest_plot[i, 6] <- "APOE4*CDK5"
}
p <- ggplot(data = forest_plot, aes(x = SF, y = Centre, xmin = lower, xmax = upper))
p <- p + geom_pointrange(aes(col = Centre))
p <- p + geom_vline(aes(fill = Centre), xintercept = 1, linetype = "dotted")
p <- p + geom_errorbar(aes(xmin = lower, xmax = upper, col = Centre),width = 0.5, cex = 1)
p <- p + facet_wrap(~gene, strip.position = "left", nrow = 1, scales = "free_y")
p <- p + geom_text(aes(label = round(p_value, 3)),hjust = -0.5, vjust = -1)
p <- p+theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank())
p <- p + coord_cartesian(xlim = c(-1, 50)) 
p




VDR_DBH_interactions <- vector(mode = "list", length = length(levels(All$Centre)))
for (centre in str_to_title(levels(All$Centre))) {
  if (centre == "Oviedo") {
    next
  }
  snp1_model <- "rs731236a"
  snp2_model <- "rs1611131r"
  covariates <- c("Sex", "Age75", "E4status")
  interaction <- paste0(snp1_model, "*", snp2_model)
  predictors <- c(covariates, interaction)
  formula <- as.formula(paste("Diag", paste(predictors, collapse = " + "), sep = " ~ "))
  args <-  list(formula = formula, family = quasibinomial("logit"), data = as.name(centre))
  linear_model <- do.call(glm, args)
  OR_summary <- Coeff.OR2(linear_model)
  cat("\n")
  print(centre)
  print(OR_summary)
  i <- which(str_to_title(levels(All$Centre)) == centre)
  VDR_DBH_interactions[[i]] <- OR_summary
  names(VDR_DBH_interactions)[i] <- centre
}

VDR_DBH_interactions <- rlist::list.clean(VDR_DBH_interactions)

#plot
forest_plot <- as.data.frame(matrix(data = NA, nrow = length(VDR_DBH_interactions), ncol = 6))
colnames(forest_plot) <- c("SF", "lower", "upper", "p_value", "Centre", "gene")
for (i in seq_along(VDR_DBH_interactions)) {
  results <- VDR_DBH_interactions[[i]]
  results <- tail(results, n=1)
  forest_plot[i, 1] <- results[3]
  forest_plot[i, 2] <- results[4]
  forest_plot[i, 3] <- results[5]
  forest_plot[i, 4] <- round(results[6], digits = 2)
  forest_plot[i, 5] <- names(VDR_DBH_interactions)[i]
  forest_plot[i, 6] <- "VDR*DBH"
}
p <- ggplot(data = forest_plot, aes(x = SF, y = Centre, xmin = lower, xmax = upper))
p <- p + geom_pointrange(aes(col = Centre))
p <- p + geom_vline(aes(fill = Centre), xintercept = 1, linetype = "dotted")
p <- p + geom_errorbar(aes(xmin = lower, xmax = upper, col = Centre),width = 0.5, cex = 1)
p <- p + facet_wrap(~gene, strip.position = "left", nrow = 1, scales = "free_y")
p <- p + geom_text(aes(label = round(p_value, 3)),hjust = -0.5, vjust = -1)
p <- p+theme(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank())


p
