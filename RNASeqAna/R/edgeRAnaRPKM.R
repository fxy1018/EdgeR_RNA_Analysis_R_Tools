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
  library(data.table)
  
  if (spe =="rat"){
    nrow = 32883  
  }
  else if (spe == "human"){
    nrow = 58233
  }
  else if (spe == "mouse"){
    nrow = 0
  }
  
  sampleDGE <- readDGE(files, path=dir, columns = c(1,3), group=group,header=FALSE, nrow = nrow)
  
  #get the gene length
  gene_length = getGeneLengths(spe)
  
  sampleDGE$genes$Length <- gene_length$length
  
  #calcculate reads per kilobase per million (RPKM) values
  # out = rpkm(sampleDGE, gene.length= gene_length$length, log=F, prior.count = 0.25)
  
  # my_rpkm = function(dge){
  #   data = dge$counts*1000*10^6
  #   print(dim(data))
  #   size = dge$samples$lib.size
  #   print(dim(size))
  #   out = apply(data, 1, function(x) x/size)
  #   out = t(out)
  #   out = apply(out, 2, function(x) x/dge$genes$Length)
  #   print(class(out))
  #   out
  # }
  
  out= rpkm(sampleDGE,normalized.lib.sizes=TRUE, prior.count = 1)
  
  # out = my_rpkm(sampleDGE)
  #write to csv file 
  outfile = file.path(dir, paste(outprefix,"_gene_expression_rpkm.csv", sep=""))
  
  write.csv(out, outfile, row.names = T, quote = F)
  
  return(out)
  
}




