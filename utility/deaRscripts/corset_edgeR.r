#!/usr/bin/env Rscript
################################################################################
### R script to compare two different conditions with count files (generated using corset and salmon) and DESeq2 packages
### Aditya Narayan Sarangi
### Designed to be executed with bulkRNASeqPIPE
################################################################################

rm(list=ls())                                        # remove all the objects from the R session
suppressMessages(library(rnaseqdea))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(DT))
suppressMessages(library(ggplot2))
suppressMessages(library(gplots))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(tximport))
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

# options list with associated default value.
option_list <- list( 
make_option(c("-P", "--projectName"),
			default=basename(getwd()),
			dest="projectName",
			help="name of the project used for storing images and tables [default: name of the current directory]."),

make_option(c("-R", "--reportName"),
			default="Corset_edgeR_HTML_Report",
			dest="reportName",
			help="name of the project used for the report [default: name of the current directory]."),


make_option(c("-t", "--targetFile"),
			default="target.txt",
			dest="targetFile",
			help="path to the design/target file [default: %default]."),

make_option(c("-T", "--templateFile"),
			dest="templateFile",
			help="path to the R markdown Template file"),


make_option(c("-q", "--corsetQuant"),
			dest="corsetQuant",
			help="path to the corset Quantification file"),

	
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
			help="median or shorth to estimate the size factors [default: %default]"),

make_option(c("-m", "--cpmCutoff"),
			default=1,
			dest="cpmCutoff", 
			help="counts-per-million cut-off to filter low counts"),
			
make_option(c("-g", "--gene.selection"),
			default="pairwise",
			dest="gene.selection", 
			help="selection of the features in MDSPlot [default: %default]"),
			
make_option(c("-n", "--normalizationMethod"),
			default="TMM",
			dest="normalizationMethod", 
			help="normalization method in calcNormFactors: \"TMM\", \"RLE\" or \"upperquartile\" [default: %default]")

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

if ( is.null(opt$corsetQuant) ) {
  stop("--path to the Corset Countfile file must be provided. See script usage (--help)")
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
corsetQuant <- opt$corsetQuant                             # path to the directory containing salmon quantification files
varInt <- opt$varInt                                 # factor of interest
condRef <- opt$condRef                               # reference biological condition
batch <- opt$batch                                   # blocking factor: NULL (default) or "batch" for example
fitType <- opt$fitType                               # mean-variance relationship: "parametric" (default), "local" or "mean"
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
locfunc <- opt$locfunc  
gene.selection <- opt$gene.selection                 # selection of the features in MDSPlot
normalizationMethod <- opt$normalizationMethod       # normalization method in calcNormFactors
cpmCutoff <- opt$cpmCutoff                             # "median" (default) or "shorth" to estimate the size factors
				

 print(paste("workDir", workDir))
 print(paste("projectName", projectName))
 print(paste("reportName", reportName))
 print(paste("targetFile", targetFile))
 print(paste("corsetQuant", corsetQuant))
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

#countData <- read.delim(corsetQuant,row.names=1)
countData <- read.delim(coresetQuant,header=TRUE,check.names=FALSE,sep="\t",row.names=1)
head (countData)

all(rownames(target) %in% colnames(countData))

countData <- countData[, rownames(target)]
all(rownames(target) == colnames(countData))

out.edgeR <- run.edgeR(counts=countData, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

dge <- out.edgeR$dge
res <- out.edgeR$results
lrt <- out.edgeR$lrt


dds = as.DESeqDataSet(dge)
#res <- results(dds)

coldata <- colData(dds)

samples <- as.factor(row.names(coldata))
samples


colData(dds) <- cbind(colData(dds), samples)

coldata <- colData(dds)
intgroup <- colnames(coldata[c(1)])
intgroup

#summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=countdata, alpha=alpha)
exportResults.edgeR(out.edgeR, group=target[,varInt], counts=countData, alpha=alpha, export=TRUE)

save.image(file=paste0(reportName, ".RData"))


report <- edgeReport(dge, lrt, projectName, intgroup, outdir = reportName,  output = 'index', nBest = 50000, nBestFeatures = 20, template = templateFile)


if(interactive()) {
    browseURL(report)
}


