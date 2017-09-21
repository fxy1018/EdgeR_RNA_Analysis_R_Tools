#' A edgeRAnaFit Function
#'
#' This function generate the fit model based on edgeR
#' @param files a dataframe which contain edgeR pathway analysis results
#' @param dir directory which count files are located
#' @param group experiment condition, which is a factor
#' @param outprefix the prefix of output file
#' @keywords pathway analysis
#' @export
#' @examples edgeRAnaFit(pathway, up=FALSE)
#' edgeRAnaFit()

edgeRAnaFit <- function(files, dir, group, outprefix, spe){
  library(edgeR)
  sampleDGE <- readDGE(files, path=dir, columns = c(1,2), group=group,header=FALSE)
  #filter the meaningful count row
  keep <- rowSums(cpm(sampleDGE) > 1) >=2
  sample_keep <- sampleDGE[keep, , keep.lib.sizes=FALSE]
  sample_keep <- calcNormFactors(sample_keep)
  #change ensemble id to entrez id

  sample_count <- as.data.frame(sample_keep$count)
  sample_count$ensembleId <- rownames(sample_keep)

  entrez_id_des <- ensemble2entrez(sample_count$ensembleId, spe)
  id_merged <- merge(sample_count, entrez_id_des, by.x="ensembleId", by.y="ensembl_gene_id", all.x=T)
  entrez_ids <- id_merged$entrezgene
  rownames(sample_keep) <- entrez_ids
  keep2 <- !is.na(entrez_ids)
  sample_keep <- sample_keep[keep2, , keep.lib.sizes=FALSE]
  sample_keep <- calcNormFactors(sample_keep)
  #Data exploration
  mds_out <- file.path(dir, paste(outprefix, "_mds.pdf", sep = ""))
  pdf(mds_out, width = 15, height =10)
  plotMDS(sample_keep)
  dev.off()

  #design matrix
  design <- model.matrix(~group)
  rownames(design) = colnames(sample_keep)

  #estimate the dispersion
  sample_keep <- estimateDisp(sample_keep, design, robust=T)

  #BCV plot
  bcv_out <- file.path(dir, paste(outprefix, "_bcv.pdf", sep = ""))
  pdf(bcv_out, width = 15, height =10)
  plotBCV(sample_keep)
  dev.off()

  fit <- glmFit(sample_keep, design)
  return(fit)

}





