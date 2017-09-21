#' A ensemble2geneLength Function
#'
#' This function allows you to convert ensemble id to entrez id based on biomart
#' @param geneIdS a vector which contain ensembl2entrez ids
#' @keywords id convert
#' @export
#' @examples ensemble2geneLength(geneIds)
#' ensemble2geneLength()



ensemble2geneLength <- function(geneIds, spe){
  
  library(biomaRt)
  mart <- useMart("ensembl")
  listDatasets(mart)
  listMarts(mart, host="www.ensembl.org")
  
  if (spe == "human"){
    mart <- useDataset("hsapiens_gene_ensembl", mart=mart)
  }else if (spe == "rat"){
    mart <- useDataset("rnorvegicus_gene_ensembl", mart=mart)
  }
  
  # geneIds = "ENSG00000000003"
  entrez_id_des <- getBM(filters= "ensembl_gene_id",
                         attributes = c("ensembl_gene_id", "start_position","end_position"),
                         values = geneIds, mart=mart)
  entrez_id_des$length = as.numeric(entrez_id_des$end_position) - as.numeric(entrez_id_des$start_position) + 1
  
  entrez_id_des
  
  return(entrez_id_des$length)
  
}

