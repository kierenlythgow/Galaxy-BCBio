
library(optparse)

#Parameters passed from bcbio.py script
option_list = list(
	make_option(c("-o", "--out_rds"), type="character", default="out.rds", 
              help="Rdata object output name [default= %default]", metavar="character"),
	make_option(c("-b", "--genome"), type="character", default=NULL,
              help="Genome build name (e.g GRCh37)", metavar="character"),
	make_option(c("-d", "--upload_dir"), type="character", default=NULL,
              help="Bcbio upload final directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genome)){
  print_help(opt_parser)
  stop("A genome build must be specified, e.g GRCh37 (genome).\n", call.=FALSE)
}

if (is.null(opt$upload_dir)){
  print_help(opt_parser)
  stop("A bcbio final directory must be supplied (upload_dir).\n", call.=FALSE)
}

#source("https://bioconductor.org/biocLite.R")
library(bcbioRNASeq)
library(data.table)
library(dplyr)
#Set regex to allow for variations of human genome builds
reg <- grepl("^grch|^hg", opt$genome, ignore.case=T)

#Check genome build to identify species
if (reg == TRUE) {
	species = 'Homo sapiens'
} else if (opt$genome == 'GRCm38') {
	species = 'Mus musculus'
}

bcb <- bcbioRNASeq(uploadDir = opt$upload_dir, organism = species, genomeBuild=opt$genome)

#Save RData object
saveRDS(bcb,file=opt$out_rds)

