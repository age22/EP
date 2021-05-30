## GENERATE SORTED TABLE OF GENES ##
df <- read.csv2(file.path(getwd(), "Input", "Epistasis_data.csv"), na.strings = "00", stringsAsFactors = FALSE) #Reading csv file and storing it as dataframe
sort(table(gsub("_.*","", colnames(df[,-(1:22)]))), decreasing = T)  #GET SORTED TABLE OF GENES OF THE SNPS PRESENT IN THE DATA


## SELECT CERTAIN GENE BELONGING SNPS AS INPUT ##
df <- read.csv2(file.path(getwd(), "Input", "Epistasis_data.csv"), 
                na.strings = c("00", "", "???", "-9"), 
                stringsAsFactors = F) 
df$DHFR_rs70991108_INDEL <- NULL
silly_function <- function(dataframe) {
  all_genes <- names(sort(table(gsub("_.*","", colnames(dataframe[,-(1:22)]))), decreasing = T))
  print(all_genes)
  cat("\n\n")
  gene <- readline(prompt = "Please select one of the previous genes to add it's snps to the input file:\t")
  snps <- grep(gene, colnames(dataframe[,-(1:22)]), value = T) # Selecting the columns with snps from the dataframe
  write(snps, file = file.path(getwd(), "Input", "snps.txt"))
}
silly_function(dataframe = df)


## SELECT ALL SNPS AS INPUT ##
df <- read.csv2(file.path(getwd(), "Input", "Epistasis_data.csv"), 
                na.strings = c("00", "", "???", "-9"), 
                stringsAsFactors = F) 
df$DHFR_rs70991108_INDEL <- NULL
for (gene in names(sort(table(gsub("_.*","", colnames(df[,-(1:22)]))), decreasing = T))) {
  snps <- grep(gene, colnames(df[,-(1:22)]), value = T)
  write(snps, file = file.path(getwd(), "Input", "snps.txt"), append = T)
}

## EXTRACT SNP FROM VCF FILES ## 

bash_script <- 'for chromosome in {1..22}
                do
                  vcftools --gzvcf \
                  ALL.chr"$chromosome".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz \
                  --snps my_snps_"$chromosome" --recode --stdout \ 
                  | grep "^"$chromosome"" | cut -d $"\t" -f 1-5 > snp_"$chromosome".txt
                done'

system(command = bash_script)

