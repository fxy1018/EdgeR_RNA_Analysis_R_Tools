#' A edgeRAnaPathview Function
#'
#' This function allows you to convert entrez id to genename based on biomart
#' @param mapped_genes a vector which contain entrez ids
#' @keywords id convert
#' @export
#' @examples
#' edgeRAnaPathview(geneIds)

edgeRAnaPathview <- function(mapped_genes, dir, spe){
  library(pathview)
  data <- mapped_genes
  pathviewFun <- function(exp, id, out, spe){
    if (spe=="human"){
      kegg_species = "hsa"
    }else if (spe=="rat"){
      kegg_species = "rno"
    }
    
    pathview(gene.data=exp, pathway.id=id, species = kegg_species, out.suffix = out)
  }
  
  conditions <- factor(data$Comparison)
  print(levels(conditions))
  if (spe=="human"){
    sub_pattern = "path:hsa"
  }else if (spe=="rat"){
    sub_pattern = "path:rno"
  }

  for (l in levels(conditions)) {
    ifelse(!dir.exists(file.path(dir, l)), dir.create(file.path(dir, l)), FALSE)
    setwd(file.path(dir, l))
    sub_data <- data[data$Comparison==l,]
    pathways <- unique(sub_data$PathID)
    for (p_id in pathways){
      path_data <- sub_data[sub_data$PathID==p_id, c("GeneID", "logFC")]
      gene_ids = path_data$GeneID
      path_data = path_data$logFC
      names(path_data) <- gene_ids
      id <- sub(sub_pattern, "", p_id)
      
      # if (id !="04723"){
      #   pathviewFun(path_data, id, l, spe)
      # }
      pathviewFun(path_data, id, l, spe)

    }

  }
  
  
}





