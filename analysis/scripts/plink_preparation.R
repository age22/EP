# Preparing input for plink

cat("Preparing input for plink...\n\n")

library(biomaRt)

## Preparing .ped file

### Transforming all factor variables into the factor datatype
factors <- -c(1,4)
All$E4status[All$E4status == 2] <- NA
All[,factors] <- lapply(All[,factors], factor)

### Getting twice the SNP calls variables, one for each allele
ped <- All[,c("ID", "Sex", "Diag", sort(rep(vector_of_snps, 2)))]

###Renaming Sex and Diag levels in a binaric way (1,2)
levels(ped$Sex) <- c("2", "1")
levels(ped$Diag) <- c("2", "1")

### For each SNP call if it is even it means is the first of the two copies (as the first SNP
### column is column number 4) of the same call, so we keep the first allele. If it is odd, 
### it means it is the second one, so we keep the second allele.

for (i in 4:(length(ped[,-c(1:3)]) + 3)) {
  if (i %% 2 == 0) {
    ped[,i] <- substr(ped[,i],1,1)
  } else {
    ped[,i] <- substr(ped[,i],2,2)
  }
}

### Coding family id, paternal id and maternal id as missing.
FID <- 0
PID <- 0
MID <- 0

### Adding and ordering the variables in the correct .ped format.

ped <- cbind(FID, PID, MID, ped)

ped_ID_Sex_Diag <- ped[c(1,4,2,3,5,6)]
ped_alleles <- ped[7:length(ped)]

ped <- cbind(ped_ID_Sex_Diag, ped_alleles)

levels(ped$Sex) <- c("2", "1", "0")
levels(ped$Diag) <- c("2", "1", "0")
ped[-c(2)][is.na(ped[-c(2)])] <- 0

factors <- -c(2)
ped[,factors] <- lapply(ped[,factors], factor)

## Preparing .map file 

### Retrieving from attribute information from Ensembl
snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
map <- data.frame(matrix(NA, nrow = length(vector_of_snps), ncol = 4))
list_of_attributes <- c("chromosome_code", "variant_id", "position_cM", "bp_coordinate")
colnames(map) <- list_of_attributes
rownames(map) <- sort(vector_of_snps)

for (snp in vector_of_snps) {
  biomart <- getBM(attributes = c('chr_name','chrom_start'),
                   filters = "snp_filter",
                   values = snp,
                   mart = snpmart)
  
  variant_id <- snp
  chromosome_code <- biomart[[1]][1]
  position_cM <- 0
  bp_coordinate <- biomart[[2]][1]
  
  for (attribute in list_of_attributes) {
    map[snp, attribute] <- get(attribute)[[1]]
  }
}

### Writing .ped and .map files
write.table(map, file = file.path(getwd(), "Input/plink", "test.map"), quote = F, sep = "\t",
            row.names = F, col.names = F)
write.table(ped, file = file.path(getwd(), "Input/plink", "test.ped"), quote = F, sep = "\t",
            row.names = F, col.names = F)
