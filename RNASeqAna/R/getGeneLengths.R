#' A getGeneLengths Function
#'
#' This function allows you to convert ensemble id to entrez id based on biomart
#' @param geneIdS a vector which contain ensembl2entrez ids
#' @keywords id convert
#' @export
#' @examples getGeneLengths()
#' getGeneLengths()



getGeneLengths <- function(spe){
  if (spe == "human"){
    data_file = 'ensembl_gene_length_human.csv'
  } else if (spe=="rat"){
    data_file = 'ensembl_gene_length_rat.csv'
  }
  dir = "C:/Users/xfan/Documents/R/my-library/RNASeqAna/R/"
  file = file.path(dir, data_file)
  out = read.csv(file)
  
  return(out)
  
}

