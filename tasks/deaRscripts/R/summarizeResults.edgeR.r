#' Summarize edgeR analysis
#'
#' Summarize edgeR analysis: diagnotic plots, dispersions plot, summary of the independent filtering, export results...
#'
#' @param out.edgeR the result of \code{run.edgeR()}
#' @param group factor vector of the condition from which each sample belongs
#' @param counts matrix of raw counts
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param col colors for the plots
#' @param log2FClim numeric vector containing both upper and lower y-axis limits for all the MA-plots produced (NULL by default to set them automatically)
#' @param padjlim numeric value between 0 and 1 for the adjusted p-value upper limits for all the volcano plots produced (NULL by default to set them automatically)
#' @return A list containing: (i) a list of \code{data.frames} from \code{exportResults.edgeR()} and (ii) a table summarizing the number of differentially expressed features
#' @author Hugo Varet

summarizeResults.edgeR <- function(out.edgeR, group, counts, alpha=0.05,
                                   col=c("lightblue","orange","MediumVioletRed","SpringGreen"),
                                   log2FClim=NULL, padjlim=NULL){  

  if (!I("tables" %in% dir())) dir.create("tables", showWarnings=FALSE)
  
  complete <- exportResults.edgeR(out.edgeR=out.edgeR, group=group, counts=counts, alpha=alpha)
  
  # small table with number of differentially expressed features
  nDiffTotal <- nDiffTotal(complete=complete, alpha=alpha)
  cat("Number of features down/up and total:\n")
  print(nDiffTotal, quote=FALSE)
  
  return(list(complete=complete))
}
