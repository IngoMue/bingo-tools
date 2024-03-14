#########################################
#                                       #
#   Calculates Dxy from mafs files      #
#                                       #
#   Author: Joshua Penalba              #
#   Date: 22 Oct 2016                   #
#   Modified by: Ingo A. Mueller        #
#   Modified Date: 12 Mar 2024          #
#                                       #
#########################################


# NOTES
# * Prior to calculating Dxy the following steps are recommended:
#   1. Run ANGSD with all populations with a -SNP_pval and -skipTriallelic flags.
#   2. Rerun ANGSD per population 
#       Use the -sites flag with a file corresponding to the recovered SNPs.
#       This will guarantee that sites with an allele fixed in one population is still included.
#       Remove the -SNP_pval flag.
#       IMPORTANT: Include an outgroup reference to polarize alleles.
#   3. Gunzip the resulting mafs files.
# 
# * Make sure the totLen only includes the chromosomes being analyzed.
# * minInd flag not added, assuming already considered in the ANGSD run.
# * Test for matching major and minor alleles not included as it would filter out sequencing errors. 
#   This has been accounted for in the allele frequency calculations.
#   This filter may give an underestimate of dxy.
# * Per site Dxy of ~0 could be common if the alternate alleles are present in a population other than the two being included in the calculation.

### Creating an argument parser
library("optparse")

option_list = list(
  make_option(c("-p","--popPair"), type="character",default=NULL,help="path to file containing desired pairwise comparisons (tab-separated)",metavar="character"),
  make_option(c("-i","--inPATH"), type="character",default=NULL,help="path in which all maf files are stored",metavar="character"),
  make_option(c("-t","--totLen"), type="numeric",default=NULL,help="total sequence length for global per site Dxy estimate [optional]",metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

### Troubleshooting input
if(is.null(opt$popPair)){
  print_help(opt_parser)
  stop("Path to input file is missing", call.=FALSE)
}

if(is.null(opt$inPATH)){
  print_help(opt_parser)
  stop("Path to input directory is missing", call.=FALSE)
}

if(is.null(opt$totLen)){
  print("Total length not supplied. The output will not be a per site estimate.")
}

### Reading data in
poppairs <- read.table(opt$popPair, sep='\t', header=T)

glob_dxy <- data.frame("Pair" = paste0(poppairs$PopA, "-", poppairs$PopB), glob_dxy = NA)

if(!is.null(opt$totLen)){
  glob_dxy$persite <- NA
}

for (n in 1:nrow(poppairs)) {

  allfreqA <- read.table(paste0(opt$inPATH, poppairs$PopA[n],".mafs"),sep='\t',row.names=NULL, header=T)
  allfreqB <- read.table(paste0(opt$inPATH, poppairs$PopB[n],".mafs"),sep='\t',row.names=NULL, header=T)
  
  ### Manipulating the table and print dxy table
  allfreq <- merge(allfreqA, allfreqB, by=c("chromo","position"))
  allfreq <- allfreq[order(allfreq$chromo, allfreq$position),]
  # -> Actual dxy calculation
  allfreq <- transform(allfreq, dxy=(knownEM.x*(1-knownEM.y))+(knownEM.y*(1-knownEM.x)))
  write.table(allfreq[,c("chromo","position","dxy")], file=paste0(poppairs$PopA[n], "_", poppairs$PopB[n], "_Dxy_persite.txt"),quote=FALSE, row.names=FALSE, sep='\t')
  cat('Created ', paste0(poppairs$PopA[n], "_", poppairs$PopB[n], "_Dxy_persite.txt\n"))
  
  ### Print global dxy
  glob_dxy$glob_dxy[n] <- sum(allfreq$dxy)
  if(!is.null(opt$totLen)){
    glob_dxy$persite[n] <- sum(allfreq$dxy)/opt$totLen
  }
}

write.table(glob_dxy, file="Dxy_global.txt",quote=FALSE, row.names=FALSE, sep='\t')
cat('Created ', "Dxy_global.txt\n")