#' A edgeRAnaRPKM Function
#'
#' This function generate the fit model based on edgeR
#' @param files a dataframe which contain edgeR pathway analysis results
#' @param dir directory which count files are located
#' @param group experiment condition, which is a factor
#' @param outprefix the prefix of output file
#' @keywords pathway analysis
#' @export
#' @examples edgeRAnaRPKM(pathway, up=FALSE)
#' edgeRAnaRPKM()

edgeRAnaRPKM <- function(files, dir, group, outprefix,spe){
  library(edgeR)
  sampleDGE <- readDGE(files, path=dir, columns = c(1,2), group=group,header=FALSE)
  
  #get the gene length
  gene_length = getGeneLengths(spe)
  
  #calcculate reads per kilobase per million (RPKM) values
  out = rpkm(sampleDGE, gene.length= gene_length$length, normalized.lib.sizes = T, log=T, prior.count = 0.25)
  
  
  #write to csv file 
  outfile = file.path(dir, paste(outprefix,"_gene_expression_log_rpkm.csv", sep=""))
  
  write.csv(out, outfile, row.names = T, quote = F)
  
  return(out)
  
}




