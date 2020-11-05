#!/usr/bin/env Rscript
################################################################################
### R script to check sanity of target file
### Adapted from SARTools 1.6.7
################################################################################
suppressMessages(library(optparse))
rm(list=ls())                                     
# options list with associated default value.
option_list <- list( 

make_option(c("-t", "--targetFile"),
			dest="targetFile",
			help="path to the design/target file"),

make_option(c("-v", "--varInt"),
			default="group",
			dest="varInt", 
			help="factor of interest [default: %default]"),

make_option(c("-c", "--condRef"),
			default="WT",
			dest="condRef",
			help="reference biological condition [default: %default]"),

make_option(c("-b", "--batch"),
			default=NULL,
			dest="batch",
			help="blocking factor [default: %default] or \"batch\" for example")
)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Program to check sanity of target file.",
					   epilogue="Computational Genome Biology Lab, IICB, Kolkata")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

# get options and arguments
workDir <- getwd()
targetFile <-opt$targetFile                       # path to the design/target file
varInt <- opt$varInt                              # factor of interest
condRef <- opt$condRef                            # reference biological condition
batch <- opt$batch                                # blocking factor: NULL (default) or "batch" for example


#############################
#Function to check sanity
#############################

checkTargetFile <- function(targetFile, varInt, condRef, batch){
  target <- read.table(targetFile, header=TRUE, sep="\t", na.strings="")
  if (!I(varInt %in% names(target))) stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target))) stop(paste("The batch effect", batch, "is not in the target file")) 
  target[,varInt] <- as.factor(target[,varInt])
  if (!I(condRef %in% as.character(target[,varInt]))) stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
  target[,varInt] <- relevel(target[,varInt],ref=condRef)
  target <- target[order(target[,varInt]),]
  rownames(target) <- as.character(target[,1])
  # check if varInt contains replicates
  if (min(table(target[,varInt]))<2) stop(paste("The factor of interest", varInt, "has a level without replicates"))
  # check if NA in the target
  if (any(is.na(cbind(target[,c(varInt, batch)], target[,1:2])))) stop("NA are present in the target file")
  # warning message if batch is numeric
  if (!is.null(batch) && is.numeric(target[,batch])) warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  if (any(grepl("[[:punct:]]", as.character(target[,varInt])))) stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
  cat("Target file:\n")
  print(target)
  return(target)
}



target <- checkTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
