#' A ensemble2entrez Function
#'
#' This function allows you to convert gene id to entrez id based on biomart
#' @param genes a vector which contain genel2entrez ids
#' @param spe sepecies
#' @keywords id convert
#' @export
#' @examples gene2entrez(geneIds)
#' gene2entrez()



ensemble2entrez <- function(genes,spe){
  library(biomaRt)
  mart <- useMart("ensembl")
  listDatasets(mart)
  listMarts(mart, host="www.ensembl.org")
  if (spe == "human"){
    mart <- useDataset("hsapiens_gene_ensembl", mart=mart)
  }
  else if (spe == "rat"){
    mart <- useDataset("rnorvegicus_gene_ensembl", mart=mart)
  }
  else if (spe == "mouse"){
    mart <- useDataset("mmusculus_gene_ensembl", mart=mart)
  }
  
  
  entrez_id_des <- getBM(filters= "external_gene_name",
                         attributes = c("ucsc", "external_gene_name", "ensembl_gene_id", "entrezgene"),
                         values = genes, mart=mart)
  
  
  entrez_id_des = entrez_id_des[duplicated(entrez_id_des$ensembl_gene_id) == F,]
  
  # atts = listAttributes(mart)
  return(entrez_id_des)
  
}



