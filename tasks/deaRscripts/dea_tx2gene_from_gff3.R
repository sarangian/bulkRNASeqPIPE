#!/usr/bin/env Rscript
suppressMessages(library(readr))
suppressMessages(library(optparse))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures))
option_list = list(
  make_option(c("-g", "--gff3"), type="character", default=NULL,help="gff3 file path", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,help="tx2gene file name (.csv)", metavar="character")
 );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
txdb <- makeTxDbFromGFF(file=opt$gff3)
txdf <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
txdf$GENEID <-paste(txdf$TXNAME,txdf$GENEID, sep = "*")
write.table(txdf,opt$out, sep=",", col.names = TRUE, row.names = FALSE, quote=FALSE)
