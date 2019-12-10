#!/usr/bin/env Rscript
suppressMessages(library(rnaseqdea))
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(plyranges))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
rm(list=ls())                                        # remove all the objects from the R session

option_list <- list( 

make_option(c("-a", "--annoType"),dest="annotation",help="annotation file type gff3 or gtf"),
make_option(c("-p", "--path"), dest="filePath",type="character", default=NULL,help="annotation file file path", metavar="character"),
make_option(c("-o", "--out"), type="character", default=NULL,help="tx2gene file name (.csv)", metavar="character")
)
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Generate tx2gene from gff3 file.",
					   epilogue="Computational Genome Biology Lab, CSIR-IICB, Kolkata 700032")    

 
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

fileType <- opt$annotation
filePath <- opt$filePath
fileOut <- opt$out   


if(fileType=="gff3"){
    tx2gene <- tr2g_gff3(filePath)
    tx2gene$gene <-paste(tx2gene$gene,tx2gene$gene_name, sep = ":")
    tx2g <-paste(tx2gene$transcript,tx2gene$gene,sep=",")
    tx2g <- str_remove_all(tx2g, "rna-")
    tx2g <- str_remove_all(tx2g, "gene-")
    colnames(tx2g) <- NULL
    names(tx2g)[1] <- "TRANSCRIPT"
    names(tx2g)[2] <- "GENE"
    write.table(tx2g,opt$out, sep=",", col.names = TRUE, row.names = FALSE, quote=FALSE)
    }

if(fileType=="gtf"){
    tx2gene <- tr2g_gtf(filePath)
    colnames(tx2gene) <- NULL    
    names(tx2gene)[1] <- "TRANSCRIPT"
    names(tx2gene)[2] <- "GENE"
    write.table(tx2gene,opt$out, sep=",", col.names = TRUE, row.names = FALSE, quote=FALSE)
    }
