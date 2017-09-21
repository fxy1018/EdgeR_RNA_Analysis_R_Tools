#' A edgeRAnaWrite Function
#'
#' This function write df Genes in txt files and return these genes
#' @param dfGenes a list of dfGene tables got from edgeRdfGenes
#' @param names a vector stroe the comparison name
#' @param dir directory which count files are located
#' @param group experiment condition, which is a factor
#' @param outprefix the prefix of output file
#' @keywords pathway analysis
#' @export
#' @examples edgeRAnaWrite()
#' edgeRAnaWrite()

edgeRAnaWrite <- function(dfGenes, comparisons, dir, outprefix){
  #get all genes' FC
  all = data.frame(EntrezID=integer(),logFC=double(), logCPM=double(),
                   LR=double(), PValue=double(),
                   FDR=double(), geneName=character() )

  for (i in (1:length(dfGenes))){
    r = dfGenes[[i]]
    data = r$table
    print(length(data))
    if (length(data) == 0){
      next
    }
    data$Comparison = rep(comparisons[[i]], dim(data)[1])
    #data$PValue <- format(data$PValue, scientific = T)
    #data$FDR <- format(data$FDR, scientific = T)
    #data <- roundDataFrame(data,3)
    data$EntrezID = rownames(data)
    data <- data[, c(8,6,7,1:5)]
    all = rbind(all,data)
    data_file <- paste(comparisons[[i]], "_fc_all.csv", sep = "")
    write.csv(data, file.path(dir,data_file), quote=F, row.names = F)
  }

  all_out = file.path(dir, paste(outprefix,"_all_dfgenes.csv", sep=""))
  write.csv(all, all_out, quote = F, row.names=F)
  return(all)
}



