#' A edgeRAnaTPM Function
#'
#' This function generate the fit model based on edgeR
#' @param files a dataframe which contain edgeR pathway analysis results
#' @param dir directory which count files are located
#' @param group experiment condition, which is a factor
#' @param outprefix the prefix of output file
#' @keywords pathway analysis
#' @export
#' @examples edgeRAnaTPM(pathway, up=FALSE)
#' edgeRAnaTPM()

edgeRAnaTPM <- function(files, dir, group, outprefix,spe){
  library(edgeR)
  sampleDGE <- readDGE(files, header=FALSE,path=dir, columns = c(1,2), group=group)
  counts = sampleDGE$counts
  #get the gene length
  gene_length = getGeneLengths(spe)
  
  lengths = gene_length$length
  count = counts[,1]
  #calcculate transcript per million (TPM) values
  tpm <- function(count, lengths) {
    rate <- (count+1) / lengths
    rate = rate / sum(rate, na.rm = T) * 1e6
    return(rate)
  }
 
  tpms <- apply(counts, 2, function(x){tpm(x, gene_length$length)})
  tpms <-as.data.frame(tpms)
  rownames(tpms)<- rownames(counts)

  #write to csv file 
  outfile = file.path(dir, paste(outprefix,"_gene_expression_tpm.csv", sep=""))
  
  write.csv(tpms, outfile, row.names = T, quote = F)
  
  return(tpms)
  
}



