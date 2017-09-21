#' A getPathway Function
#'
#' This function allows you to get pathwasy with q-value <0.05
#' @param pathways a dataframe which contain edgeR pathway analysis results
#' @param fdr cutoff of pathway statistic
#' @param out whether to write the result to files
#' @param dir output directory
#' @param outprefix prefix of output files
#' @keywords pathway analysis
#' @export
#' @examples getPathway(pathway, up=FALSE)
#' getpathway()

getPathway <- function(pathways, fdr=0.05, out = T, dir = NA, outprefix =NA){
  #get significatn pathways
  up <-pathways[pathways$FDR_up < fdr,]
  up$regulation <- rep("UP", dim(up)[1])
  down <- pathways[pathways$FDR_down < fdr,]
  down$regulation <- rep("DOWN", dim(down)[1])

  sig_path <- rbind(up, down, make.row.names=F)
  if (out==T){
    out_file = paste(outprefix, "_kegg_", toString(fdr),".txt", sep="")
    print(paste("writing", out_file, "......"))
    write.table(sig_path, file.path(dir,out_file), sep="\t", quote=F)
  }
  return(sig_path)
}

