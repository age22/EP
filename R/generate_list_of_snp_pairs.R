#' Generates a list of snp pairs to be used in further analysis steps
#'
#' Generates a list of snp pairs with (reported_model) or without (best_model)
#' the snp model.
#'
#' @param mode string, whether it should use a snp model in the interactions
#' list ("reported_model") or just the snp in question ("best_model")
#' @param snp_df dataframe, the snp with the snp interactions with the variable
#' "snp1" being the first member of the interaction and "snp2" the second. If
#' using models, "model1" should define the first member of the interaction and
#' model2 the second one.
#'
#' @return A list of character vectors, containing the character vectors the snp
#' names or snp models interactions, depending on the mode of operation selected
#' @export
generate_list_of_snp_pairs <-
  function(mode, snp_df) {

    nrows <- nrow(snp_df)
    list_of_snp_pairs <- vector(mode = "list", length = nrows)

    if (mode == "best_model") {
      for (i in seq_len(nrows)) {
        snp_1 <- snp_df$snp1[i]
        snp_2 <- snp_df$snp2[i]
        snp_pair <- c(snp_1, snp_2)
        list_of_snp_pairs[[i]] <- snp_pair
      }
    }

    if (mode == "reported_model") {
      for (i in seq_len(nrows)) {
        snp_1 <- snp_df$model1[i]
        snp_2 <- snp_df$model2[i]
        snp_pair <- c(snp_1, snp_2)
        list_of_snp_pairs[[i]] <- snp_pair
      }
    }
    list_of_snp_pairs <- unique(list_of_snp_pairs)
  }
