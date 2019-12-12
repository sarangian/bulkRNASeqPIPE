#!/usr/bin/env Rscript
################################################################################
### R script to compare two different conditions with count files (generated using featureCounts) and DESeq2 packages
### Aditya Narayan Sarangi
### Designed to be executed with bulkRNASeqPIPE
################################################################################

rm(list=ls()) 
suppressMessages(library(rnaseqdea))
suppressMessages(library(DESeq2))
suppressMessages(library(DT))
suppressMessages(library(ggplot2))
suppressMessages(library(gplots))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(regionReport))
suppressMessages(library(DEFormats))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(dplyr))
suppressMessages(library(colorspace))
suppressMessages(library(optparse))
suppressMessages(library(scales))
suppressMessages(library(readr))
suppressMessages(library(DelayedArray))
                                 # to run the script in command lines

# options list with associated default value.
option_list <- list( 
make_option(c("-P", "--projectName"),
			default=basename(getwd()),
			dest="projectName",
			help="name of the project used for storing images and tables [default: name of the current directory]."),

make_option(c("-R", "--reportName"),
			default="featureCounts_DESeq2_HTML_Report",
			dest="reportName",
			help="name of the project used for the report [default: name of the current directory]."),


make_option(c("-t", "--targetFile"),
			default="target.txt",
			dest="targetFile",
			help="path to the design/target file [default: %default]."),

make_option(c("-T", "--templateFile"),
			dest="templateFile",
			help="path to the R markdown Template file"),


make_option(c("-q", "--featureQuant"),
			dest="featureQuant",
			help="path to the featureCount Quantification file [default: %default]."),

	
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
			help="blocking factor [default: %default] or \"batch\" for example"),

make_option(c("-f", "--fitType"),
			default="parametric",
			dest="fitType", 
			help="mean-variance relationship: [default: %default],local or mean"),

make_option(c("-a", "--alpha"),
			default=0.05,
			dest="alpha", 
			help="threshold of statistical significance [default: %default]"),

make_option(c("-p", "--pAdjustMethod"),
			default="BH",
			dest="pAdjustMethod", 
			help="p-value adjustment method: \"BH\" or \"BY\" [default: %default]"),


make_option(c("-l", "--locfunc"),
			default="median",
			dest="locfunc", 
			help="median or shorth to estimate the size factors [default: %default]")

)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Compare two or more biological conditions in a RNA-Seq framework with DESeq2.",
					   epilogue="Computational Genome Biology Lab, CSIR-IICB, Kolkata 700032")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options



#Check mandetory inputs 
if ( is.null(opt$targetFile) ) {
  stop("--sample groupfile file / target file must be provided. See script usage (--help)")
}

if ( is.null(opt$featureQuant) ) {
  stop("--path to the featureCounts Countfile file must be provided. See script usage (--help)")
}


if ( is.null(opt$condRef) ) {
  stop("--reference biological condition name must be provided. See script usage (--help)")
}





# get options and arguments
workDir <- getwd()
projectName <- opt$projectName 
reportName <-opt$reportName                      # name of the project
targetFile <- opt$targetFile   
templateFile <- opt$templateFile               # path to the design/target file   
featureQuant <- opt$featureQuant                            # path to the directory containing salmon quantification files
varInt <- opt$varInt                                 # factor of interest
condRef <- opt$condRef                               # reference biological condition
batch <- opt$batch                                   # blocking factor: NULL (default) or "batch" for example
fitType <- opt$fitType                               # mean-variance relationship: "parametric" (default), "local" or "mean"
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
locfunc <- opt$locfunc                               # "median" (default) or "shorth" to estimate the size factors
				

 print(paste("workDir", workDir))
 print(paste("projectName", projectName))
 print(paste("reportName", reportName))
 print(paste("targetFile", targetFile))
 print(paste("featureQuant", featureQuant))
 print(paste("varInt", varInt))
 print(paste("condRef", condRef))
 print(paste("batch", batch))
 print(paste("fitType", fitType))
 print(paste("alpha", alpha))
 print(paste("pAdjustMethod", pAdjustMethod))
 print(paste("locfunc", locfunc))

################################################################################
###                             running script                               ###
################################################################################

dir.create("tables", showWarnings = FALSE, recursive = TRUE)

					   
# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
rownames(target) <- as.character(target[,1])

countData <- read.table(featureQuant, header=TRUE, row.names=1)
countData <- countData[ ,6:ncol(countData)]
colnames(countData) <- gsub("\\.[sb]am$", "", colnames(countData))
colnames(countData) <- lapply(colnames(countData), function(x) sapply(strsplit(x, "\\."), tail, 1))
colnames(countData)
countData <- as.matrix(countData)


all(rownames(target) %in% colnames(countData))

countData <- countData[, rownames(target)]
all(rownames(target) == colnames(countData))

#analysis with DESeq2
out.DESeq2 <- run.DESeq2_corset(counts=countData, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         alpha=alpha)
#####
dds <- out.DESeq2$dds
res <- results(dds)


dds <- dds[ rowSums(counts(dds)) > 0, ]


coldata <- colData(dds)
intgroup<- colnames(coldata[c(2,3)])


dds.rld.trans <- rlog(dds, blind=FALSE)


sampleDists <- as.matrix(dist(t(assay(dds.rld.trans))))

exportResults.DESeq2(out.DESeq2, group=unique(target$group), alpha=alpha)

save.image(file=paste0(reportName, ".RData"))

report <- DESeq2Report(dds, projectName, intgroup, outdir = reportName,  output = 'index', nBest = 50000, nBestFeatures = 20, digits = 3, template = templateFile)

if(interactive()) {
    browseURL(report)
}


