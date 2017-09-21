
#' A ensemble2entrez Function
#'
#' This function allows you to convert ensemble id to entrez id based on biomart
#' @param geneIdS a vector which contain ensembl2entrez ids
#' @keywords id convert
#' @export
#' @examples ensemble2entrez(geneIds)
#' ensemble2entrez()



ensemble2entrez <- function(geneIds,spe){
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
  

  entrez_id_des <- getBM(filters= "ensembl_gene_id",
                         attributes = c("ensembl_gene_id", "entrezgene"),
                         values = geneIds, mart=mart)
  entrez_id_des = entrez_id_des[duplicated(entrez_id_des$ensembl_gene_id) == F,]


  return(entrez_id_des)

}













