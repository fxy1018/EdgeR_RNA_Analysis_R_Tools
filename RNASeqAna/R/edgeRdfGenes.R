#' A edgeRdfGenes Function
#'
#' This function allows you to get pathwasy with q-value <0.05
#' @param fit fit object generated by edgeRAnaFit
#' @param coef pairwise group
#' @keywords pathway analysis
#' @export
#' @examples edgeRdfGenes(pathway, up=FALSE)
#' edgeRdfGenes()


edgeRdfGenes <- function(fit, coef, pvalue=1,spe){
  if (length(coef) == 1) {
    lrt <- glmLRT(fit, coef=coef)
  }
  else{
    lrt <- glmLRT(fit, contrast=coef)
  }
  #get all genes
  gene <- topTags(lrt, p.value= pvalue,n = Inf)
  gene_ids <- rownames(gene)
  if (length(gene_ids) != 0 ){
    gene_name <- entrez2genename(gene_ids, spe)
    gene_names <- c()
    #let converted gene_name has the same order as gene_ids, the biomart convertion require input is sorted ids
    #the converted names may have no name, 1 name, and multiple names, moidfied it.
    for (i in (1:length(gene_ids))){
      id <- gene_ids[i]
      name <- gene_name$external_gene_name[gene_name$entrezgene==id]

      if (length(name) == 0){
        name <- NA
      }
      else if (length(name) > 1){
        name <- name[1]
      }
      gene_names <- c(gene_names, name)
    }
    gene$table$gene <- gene_names
  }
  out <- list(lrt, gene)
  return(out)
}

