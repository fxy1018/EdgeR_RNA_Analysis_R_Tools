#' A edgeRAnaReactome Function
#'
#' This function is to did reactome pathway analysis
#' @param df a dataframe with gene differential expression genes
#' @param fdr q value cutoff for differentially expressed genes. Numeric value between 0 and 1
#' @param dir directory which count files are located
#' @param comparisons pairwise comparisons
#' @param outprefix the prefix of output file
#' @keywords pathway analysis
#' @export
#' @examples edgeRAnaReactome(pathway, up=FALSE)
#' edgeRAnaReactome()

edgeRAnaReactome <- function(df, dfCutOff= 0.05, fdr=0.99, dir, comparisons, outprefix, spe){
  library(ReactomePA)
  all = data.frame(PathwayID=character(), Description=character(),
                   GeneRatio=character(), BgRatio=character(),
                   pvalue=double(), p.adjust = double(),
                   geneID=character(), Count=integer(),
                   Comparison=character())
  if (spe=="human"){
    react_species = "human"
  } else if (spe=="rat"){
    react_species = "rat"
  }
  
  
  for (c in comparisons){
    #get significant df gene for each comparison
    tmp_df <- df[df$Comparison==c & df$FDR <dfCutOff,]
    if (dim(tmp_df)[1] == 0){
      next
    }
    
    if (is.data.frame(tmp_df)) {
      tmp_df$geneID <- lapply(tmp_df$EntrezID, toString)
      
      y <- enrichPathway(gene=tmp_df$geneID,  organism = react_species, pAdjustMethod='BH', qvalueCutoff=fdr, readable=T)
      if (!is.null(y)){
        y@result$Comparison <- rep(c, dim(y)[1])
        all = rbind(all, y@result)
        
        filename <- paste(c,toString(fdr),"reactome.txt", sep="_")
        filepath <- file.path(dir, filename)
        write.table(y, filepath, sep="\t", row.names = F, quote = F)
      }
  }}
  
 
  all_out = file.path(dir, paste(outprefix, "_reactome_all.txt", sep=""))
  write.table(all, all_out, sep="\t", quote=F, row.names = F)
  return(all)
}





