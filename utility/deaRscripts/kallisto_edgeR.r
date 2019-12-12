#!/usr/bin/env Rscript
################################################################################
### R script to compare two different conditions with kallisto and edgeR packages
### Aditya Narayan Sarangi
### Designed to be executed with bulkRNASeqPIPE
################################################################################
rm(list=ls())                      # remove all the objects from the R session
suppressMessages(library(rnaseqdea))
suppressMessages(library(edgeR))
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
suppressMessages(library(rhdf5))
suppressMessages(library(tximport))suppressMessages(library(rhdf5))
suppressMessages(library(tximport))
suppressMessages(library(DelayedArray))

# options list with associated default value.
option_list <- list( 
make_option(c("-P", "--projectName"),
			default=basename(getwd()),
			dest="projectName",
			help="name of the project used for storing images and tables [default: name of the current directory]."),

make_option(c("-R", "--reportName"),
			default="Kallisto_edgeR_HTML_Report",
			dest="reportName",
			help="name of the project used for the report [default: name of the current directory]."),

make_option(c("-T", "--templateFile"),
			dest="templateFile",
			help="path to the R markdown Template file"),

make_option(c("-t", "--targetFile"),
			default="target.txt",
			dest="targetFile",
			help="path to the design/target file [default: %default]."),

make_option(c("-q", "--quantDir"),
			dest="quantDir",
			help="path to the directory containing the Kallisto Quantification Folder"),		

make_option(c("-G", "--tx2geneDirectory"),
			default="NULL",
			dest="tx2geneDirectory",
			help="path to the tx2gene Directory [default: %default]."),
		
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

make_option(c("-a", "--alpha"),
			default=0.05,
			dest="alpha", 
			help="threshold of statistical significance [default: %default]"),

make_option(c("-p", "--pAdjustMethod"),
			default="BH",
			dest="pAdjustMethod", 
			help="p-value adjustment method: \"BH\" or \"BY\" [default: %default]"),

make_option(c("-m", "--cpmCutoff"),
			default=1,
			dest="cpmCutoff", 
			help="counts-per-million cut-off to filter low counts"),
				
make_option(c("-n", "--normalizationMethod"),
			default="TMM",
			dest="normalizationMethod", 
			help="normalization method in calcNormFactors: \"TMM\", \"RLE\" or \"upperquartile\" [default: %default]")
)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Compare two or more biological conditions in a RNA-Seq framework with edgeR.",
					   epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

# get options and arguments
workDir <- getwd()
projectName <- opt$projectName  
reportName <- opt$reportName                     # name of the project
targetFile <- opt$targetFile  
templateFile <-opt$templateFile                       # path to the design/target file
quantDir <- opt$quantDir							 # path to the directory containing raw counts files
tx2geneDirectory <- opt$tx2geneDirectory		 
varInt <- opt$varInt                                 # factor of interest
condRef <- opt$condRef                               # reference biological condition
batch <- opt$batch                                   # blocking factor: NULL (default) or "batch" for example
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
normalizationMethod <- opt$normalizationMethod       # normalization method in calcNormFactors
cpmCutoff <- opt$cpmCutoff                           # counts-per-million cut-off to filter low counts

print(paste("projectName", projectName))
print(paste("targetFile", targetFile))
print(paste("varInt", varInt))
print(paste("condRef", condRef))
print(paste("batch", batch))
print(paste("alpha", alpha))
print(paste("pAdjustMethod", pAdjustMethod))
print(paste("normalizationMethod", normalizationMethod))
print(paste("templateFile", templateFile))

################################################################################
###                             running script                               ###
################################################################################

#Check mandetory inputs 
if ( is.null(opt$targetFile) ) {
  stop("--sample groupfile file / target file must be provided. See script usage (--help)")
}

if ( is.null(opt$quantDir) ) {
  stop("--path to the directory containing the Kallisto Quantification files must be provided. See script usage (--help)")
}

if ( is.null(opt$tx2geneDirectory) ) {
  stop("--tx2geneDirectory folder path must be provided. See script usage (--help)")
}
if ( is.null(opt$condRef) ) {
  stop("--reference biological condition name must be provided. See script usage (--help)")
}

################################################################################

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
				   
# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
#names(target$files)
print (target)
group=unique(target[,varInt])
group

#######################################
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)
#names(target$files)
print (target)
sampledata = read.table(targetFile, header = TRUE, sep='\t')
print (sampledata)
samplenames = sampledata$samples
print (samplenames)
files <- file.path(quantDir, 'kallisto_quant', target$samples, "abundance.h5")

kallisto_dname <- dirname(files)
kallisto_file_names <- basename(kallisto_dname)
names(files) <- kallisto_file_names
print (files)
tx2gene <- read_csv(file.path(tx2geneDirectory, "tx2gene.csv"))
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene)
counts <- txi.kallisto$counts 
head(counts) 



#######################################
rownames_countMat <- rownames(counts)

#print (rownames_countMat)

#str(counts)
#countdata = as.data.frame(lapply(counts, as.integer))
#str(countdata)

write.table(counts, file=paste0(quantDir,"counts.txt"), row.names=FALSE, sep="\t", quote=FALSE)

DF <- read.table(file=paste0(quantDir,"counts.txt"), sep = "\t", header=TRUE, quote="", fill=FALSE)

head (DF)
countdata = as.data.frame(lapply(DF, as.integer))

row.names(countdata) <- rownames_countMat

str(countdata)
head(countdata)


out.edgeR <- run.edgeR_trans(counts=countdata, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)


dge <- out.edgeR$dge
res <- out.edgeR$results
lrt <- out.edgeR$lrt

dds = as.DESeqDataSet(dge)
#res <- results(dds)

coldata <- colData(dds)

samples <- as.factor(row.names(coldata))

colData(dds) <- cbind(colData(dds), samples)

coldata <- colData(dds)
intgroup <- colnames(coldata[c(1)])

exportResults.edgeR(out.edgeR, group=group, counts=countdata, alpha=alpha, export=TRUE)

save.image(file=paste0(reportName, ".RData"))

report <- edgeReport(dge, lrt, projectName, intgroup, outdir = reportName,  output = 'index', nBest = 50000, nBestFeatures = 20, template = templateFile)


if(interactive()) {
    browseURL(report)
}
