#' A entrez2genename Function
#'
#' This function allows you to convert entrez id to genename based on biomart
#' @param geneIdS a vector which contain entrez ids
#' @keywords id convert
#' @export
#' @examples
#' entrez2genename(geneIds)




entrez2genename <- function(geneIds, spe){
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


  gene_names <- getBM(filters= "entrezgene",
                         attributes = c("entrezgene", "external_gene_name","description"),
                         values = geneIds, mart=mart)
  gene_names = gene_names[duplicated(gene_names$entrezgene) == F,]

  return(gene_names)

}













