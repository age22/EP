
# CENTRES IDS
for (centre in levels(All$Centre)) {
  ids <- All$ID[All$Centre == centre]
  write(ids, file = paste(centre, "ids.txt", sep = "_"), sep = "\n")
}

#REGIONS IDS
for (region in names(table(All$Region))) {
  ids <- All$ID[All$Region == region]
  write(ids, file = paste(region, "ids.txt", sep = "_"), sep = "\n")
}

vector_of_centres <- vector(length = length(pop_gpa$Centre))
for (i in seq_len(nrow(pop_gpa))) {
 id <- pop_gpa[i,]$IID
 id <- paste0("^", id, "$")
 x <- grep(id, All$ID)
 centre <- as.character(All[x,]$Centre)
 vector_of_centres[i] <- centre
}
