
library(optparse)

#Parameters for bioc_coldata
option_list = list(
	make_option(c("-g", "--ignore"), action="store_false", default=NULL,
              help="ignore mismatches"),
	make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="metadata file name", metavar="character"),
	make_option(c("-i", "--in_rds"), type="character", default=NULL,
              help="Pre-saved Rdata object input name", metavar="character"),
	make_option(c("-o", "--joined_rds"), type="character", default="joined.rds", 
              help="Joined Rdata object output name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Check parameters for metadata file and RData object are specified
if (is.null(opt$metadata)){
  print_help(opt_parser)
  stop("A metadata file must be supplied.\n", call.=FALSE)
}

if (is.null(opt$in_rds)){
  print_help(opt_parser)
  stop("An RData object file must be supplied.\n", call.=FALSE)
}


#source("https://bioconductor.org/biocLite.R")
library(bcbioRNASeq)
library(data.table)
library(dplyr)

# Import data and combine columns based on rownames (sample names)
bcb <- readRDS(opt$in_rds)
pheno <- data.frame(colData(bcb))
meta <- read.table(opt$metadata, header=TRUE)
	
if (!is.null(opt$ignore)){
	print("Ignoring sample mismatches")

	#Merge the bcbio data with the metadata using inner joins
	mergeddf <- merge(pheno, meta, by.x = 'sampleName', by.y = 'Sample', all=FALSE)
	
	#Check that at least 2 samples are in the dataset
	if (nrow(mergeddf) < 2) {

		stop("There are less than 2 samples so differential analyis is not possible")

	} else {
		cat("\nMerging the sequencing QC metrics & phenotypic data  ...  ")
        
		#Subset the data based on the specific matched sample ids
		bcb <- bcb[,bcb$sampleName %in% meta$Sample]
		colData(bcb) <- DataFrame(mergeddf)
        	print(colData(bcb))
	}

} else {
	unmatched = anti_join(pheno, meta, by = c("sampleName" = "Sample"))
	
	#Check their are no mismatches after antijoining
	if (nrow(unmatched) == 0){
		print("No sample mismatches")
		cat("\nMerging the sequencing QC metrics & phenotypic data  ...  ")
	
		#Merge table using inner joins
		mergeddf <- merge(pheno, meta, by.x = 'sampleName', by.y = 'Sample', all=FALSE)
		colData(bcb) <- DataFrame(mergeddf)
		print(colData(bcb))
	} else {
		print(unmatched)
		stop("There are sample mismatches")
	}
}

#Save RData object
saveRDS(bcb,file=opt$joined_rds)

