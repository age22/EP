check <- unlist(Filter(Negate(is.null), sapply(list_of_objects, function(SNP){if(SNP$minor_allele != SNP$imput_minor_allele){SNP$id}})), use.names = F)
