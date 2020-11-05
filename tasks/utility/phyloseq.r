#!/usr/bin/env Rscript
################################################################################
### R script to compare two different conditions with count file (generated using corset and salmon) and DESeq2 package
### Aditya Narayan Sarangi
### Designed to be executed with bulkRNASeqPIPE
################################################################################

rm(list=ls())                                        # remove all the objects from the R session
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(gplots))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(pheatmap))
suppressMessages(library(optparse))
suppressMessages(library(phyloseq))
suppressMessages(library(ampvis2))
# to run the script in command lines

# options list with associated default value.
option_list <- list( 
make_option(c("-o", "--outFolder"),
			dest="outFolder",
			help="name of the project used for storing images and tables"),

make_option(c("-t", "--targetFile"),
			dest="targetFile",
			help="path to the mapping file [default: %default]."),

make_option(c("-P", "--phyloseqBIOM"),
			dest="phyloseqBIOM",
			help="path to the corset Quantification file [default: %default]."),

make_option(c("-A", "--ampvisOTU"),
			dest="ampvisOTU",
			help="path to the corset Quantification file [default: %default]."),
	
make_option(c("-v", "--varInt"),
			dest="varInt", 
			help="factor of interest"),

make_option(c("-c", "--condRef"),
			default="WT",
			dest="condRef",
			help="reference biological condition [default: %default]"),

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
					   epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options


#Check mandetory inputs 
if ( is.null(opt$targetFile) ) {
  stop("--sample groupfile file / target file must be provided. See script usage (--help)")
}

if ( is.null(opt$phyloseqBIOM) ) {
  stop("--path to the phyloseq biom file must be provided. See script usage (--help)")
}


# get options and arguments
workDir <- getwd()
projectName <- opt$projectName 
reportName <-opt$reportName                      	 # name of the project
targetFile <- opt$targetFile  
templateFile <-opt$templateFile                      # path to the design/target file   
ampvisOTU <- opt$ampvisOTU                     		 # path to the directory containing salmon quantification files
phyloBIOM <-opt$phyloseqBIOM
varInt <- opt$varInt                                 # factor of interest
condRef <- opt$condRef                               # reference biological condition
fitType <- opt$fitType                               # mean-variance relationship: "parametric" (default), "local" or "mean"
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
locfunc <- opt$locfunc                               # "median" (default) or "shorth" to estimate the size factors
				

print(paste("workDir", workDir))
print(paste("projectName", projectName))
print(paste("reportName", reportName))
print(paste("targetFile", targetFile))
print(paste("phyloQUANT", phyloBIOM))
print(paste("ampvisOTU", ampvisOTU))
print(paste("varInt", varInt))
print(paste("condRef", condRef))
print(paste("fitType", fitType))
print(paste("alpha", alpha))
print(paste("pAdjustMethod", pAdjustMethod))
print(paste("locfunc", locfunc))

################################################################################
###                             running script                               ###
################################################################################

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
otu<-import_biom(phyloBIOM)
map<-import_qiime_sample_data(targetFile)
sample_names(otu)
sample_names(map)
data <- merge_phyloseq(otu, map)
save(data,file="Taxa.phyloseq")
ntaxa(data)
rank_names(data)
sample_variables(data)

colnames(tax_table(data)) = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
#colnames(tax_table(data))



#dds<-DESeq(dds,test="Wald",fitType=fitType)

#intgroup
#res <-results(dds)


############

run.DESeq2 <- function(counts, targetFile, varInt,
                       locfunc="median", fitType="parametric", pAdjustMethod="BH",alpha=0.05, ...){
  # building dds object
  dds = phyloseq_to_deseq2(data, design=formula(paste("~",varInt)))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)),collapse=" "),"\n")           
  
  # normalization
  dds <- estimateSizeFactors(dds,locfunc=eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  
  # estimating dispersions
  dds<-DESeq(dds,test="Wald",fitType=fitType)  # statistical testing: perform all the comparisons between the levels of varInt

  results <- list()
  for (comp in combn(nlevels(colData(dds)[,varInt]), 2, simplify=FALSE)){
    levelRef <- levels(colData(dds)[,varInt])[comp[1]]
    levelTest <- levels(colData(dds)[,varInt])[comp[2]]
    results[[paste0(levelTest,"_vs_",levelRef)]] <- results(dds, contrast=c(varInt, levelTest, levelRef),
                                                            pAdjustMethod=pAdjustMethod, alpha=alpha)
    cat(paste("Comparison", levelTest, "vs", levelRef, "done\n"))
  }
  
  return(list(dds=dds,results=results,sf=sizeFactors(dds)))
}


out.DESeq2 <- run.DESeq2(counts=data, target=target, varInt=varInt,fitType=fitType, pAdjustMethod=pAdjustMethod,alpha=alpha)
#####
dds <- out.DESeq2$dds
res <- results(dds)

coldata<-colData(dds)
#print(coldata)
intgroup <-(unique(coldata[, 2]))
############################
#EXPORT RESULT
exportResults.DESeq2 <- function(out.DESeq2, group, alpha=0.05, export=TRUE){
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # comptages bruts et normalis?s
  counts <- data.frame(Id=rownames(counts(dds)), counts(dds), round(counts(dds, normalized=FALSE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
  # baseMean avec identifiant
  bm <- data.frame(Id=rownames(results[[1]]),baseMean=round(results[[1]][,"baseMean"],2))
  # merge des info, comptages et baseMean selon l'Id
  base <- merge(counts, bm, by="Id", all=TRUE)
  tmp <- base[,paste("norm", colnames(counts(dds)), sep=".")]
  for (cond in levels(group)){
    base[,cond] <- round(apply(as.data.frame(tmp[,group==cond]),1,mean),0)
  }
  
  complete <- list()
  for (name in names(results)){
    complete.name <- base

    # ajout d'elements depuis results
    res.name <- data.frame(Id=rownames(results[[name]]),
                           FoldChange=round(2^(results[[name]][,"log2FoldChange"]), 3),
                           log2FoldChange=round(results[[name]][,"log2FoldChange"], 3),
                           stat=round(results[[name]][,"stat"], 3),
                           pvalue=results[[name]][,"pvalue"],
                           padj=results[[name]][,"padj"])
    complete.name <- merge(complete.name, res.name, by="Id", all=TRUE)
    # ajout d'elements depuis mcols(dds)
    mcols.add <- data.frame(Id=rownames(counts(dds)),dispGeneEst=round(mcols(dds)$dispGeneEst,4),
                            dispFit=round(mcols(dds)$dispFit,4),dispMAP=round(mcols(dds)$dispMAP,4),
                            dispersion=round(mcols(dds)$dispersion,4),betaConv=mcols(dds)$betaConv,
                            maxCooks=round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by="Id", all=TRUE)
    complete[[name]] <- complete.name
    
    if (export){
      # s?lection des up et down
      up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0),]
      up.name <- up.name[order(up.name$padj),]
      down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange<=0),]
      down.name <- down.name[order(down.name$padj),]
      
      # exports
      name <- gsub("_","",name)
      write.table(complete.name, file=paste0("tables/",name,".complete.txt"), sep="\t", row.names=FALSE, dec=".", quote=FALSE)
      write.table(up.name, file=paste0("tables/", name,".up.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
      write.table(down.name, file=paste0("tables/", name,".down.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
    }
  }

  return(complete)
}


exportResults.DESeq2(out.DESeq2,group=intgroup, alpha=alpha)



############################


#Volcano plot
volcano_plot_FC_15_P_05<-EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-8, 8),
    title = 'FC Cutoff 1.5, P-value cutoff 0.05',
    pCutoff = 0.05,
    FCcutoff = 1.5,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1)

ggsave(filename = "volcano_plot_FC_15_P_05.tiff", volcano_plot_FC_15_P_05, width = 18, dpi = 300, units = "cm", device='tiff')

#Volcano plot
volcano_plot_FC_2_P_05<-EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim = c(-8, 8),
    title = 'FC Cutoff 2, P-value cutoff 0.05',
    pCutoff = 0.05,
    FCcutoff = 2,
    col=c('black', 'black', 'black', 'red3'),
    colAlpha = 1)

ggsave(file=paste0("figures/","volcano_plot_FC_2_P_05.tiff"), volcano_plot_FC_2_P_05, width = 18, dpi = 300, units = "cm", device='tiff')

res = res[order(res$padj, na.last=NA), ]
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
#head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

differential_abundance_genus <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(file = paste0("figures/","differential_abundance_genus.tiff"), differential_abundance_genus, width = 18, dpi = 300, units = "cm", device='tiff')


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
x = sort(x, TRUE)
sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))

differential_abundance_species <- ggplot(sigtab, aes(x=Species, y=log2FoldChange, color=Phylum)) + geom_point(size=6) +
theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(file = paste0("figures/","differential_abundance_species.tiff"), differential_abundance_species, width = 18, dpi = 300, units = "cm", device='tiff')

#plot ampvis2 Species
data <- amp_load(otutable = ampvisOTU, metadata = targetFile)
plot_box_species<-amp_boxplot(data,
					tax_aggregate = "Species",
					group_by = varInt,
					tax_show = 25,
					tax_add = "Phylum", 
					normalise = FALSE, 
					plot_log =TRUE, 
					point_size = 0.1, 
					sort_by = "median") + geom_boxplot(alpha =0.05,outlier.size=0.5) + scale_color_manual(values=c("red","blue"))

ggsave(file = paste0("figures/","plot_box_species.tiff"), plot_box_species, width = 18, dpi = 300, units = "cm", device='tiff')


#AMPVIS_GENUS

plot_box_genus<-amp_boxplot(data,
					tax_aggregate = "Genus",
					group_by = varInt,
					tax_show = 25,
					tax_add = "Phylum", 
					normalise = FALSE, 
					plot_log =TRUE, 
					point_size = 0.1, 
					sort_by = "median") + geom_boxplot(alpha =0.05,outlier.size=0.5) + scale_color_manual(values=c("red","blue"))

ggsave(file = paste0("figures/","box_plot_genus.tiff"), plot_box_genus, width = 18, dpi = 300, units = "cm", device='tiff')

#plot ampvis2 Family
plot_box_family<-amp_boxplot(data,
					tax_aggregate = "Family",
					group_by = varInt,
					tax_show = 25,
					tax_add = "Phylum", 
					normalise = FALSE, 
					plot_log =TRUE, 
					point_size = 0.1, 
					sort_by = "median") + geom_boxplot(alpha =0.05,outlier.size=0.5) + scale_color_manual(values=c("red","blue"))

ggsave(file = paste0("figures/","box_plot_family.tiff"), plot_box_family, width = 18, dpi = 300, units = "cm", device='tiff')


#plot species heatmap
plot_heatmap_species <- amp_heatmap(data,
            group_by = varInt,
            tax_aggregate = "Species",
            tax_add = "Phylum",
            tax_show = 25,
			normalise = FALSE,
            color_vector = c("white", "darkred"),
            plot_colorscale = "sqrt",
            plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.position="right")

ggsave(file = paste0("figures/","species_heatmap.tiff"), plot_heatmap_species, width = 18, dpi = 300, units = "cm", device='tiff')


#plot genus heatmap
plot_heatmap_genus <- amp_heatmap(data,
            group_by = varInt,
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            tax_show = 25,
			normalise = FALSE,
            color_vector = c("white", "darkred"),
            plot_colorscale = "sqrt",
            plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.position="right")
ggsave(file = paste0("figures/","genus_heatmap.tiff"), plot_heatmap_genus, width = 18, dpi = 300, units = "cm", device='tiff')

#plot genus heatmap
plot_heatmap_family <- amp_heatmap(data,
            group_by = varInt,
            tax_aggregate = "Family",
            tax_add = "Phylum",
            tax_show = 25,
			normalise = FALSE,
            color_vector = c("white", "darkred"),
            plot_colorscale = "sqrt",
            plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
        axis.text.y = element_text(size=8),
        legend.position="right")
ggsave(file = paste0("figures/","family_heatmap.tiff"), plot_heatmap_family, width = 18, dpi = 300, units = "cm", device='tiff')