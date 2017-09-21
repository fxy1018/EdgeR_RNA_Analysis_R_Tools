#' A qcPlots Function
#'
#' This function allows you to draw a pca plot
#' @param files the files came from htseq-count (only contain two coloumns)
#' @param dir directory which count files are located
#' @param condition experiment condition, which is a factor
#' @param outprefiex output file name
#' @keywords qc plot
#' @export
#' @examples qcPlot(files, dir, condtion, outprefix)
#' qcPlot()

qcPlots <- function(files, dir, condition, outprefix){
  library("vsn")
  library("ggplot2")
  library("DESeq2")
  sampleTable <- data.frame(sampleName = files,
                            fileName = files,
                            condition = condition)

  dds<- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                   directory = dir,
                                   design= ~ condition)

  #filter the data
  dds_filtered <- dds[ rowSums(counts(dds)) > 1, ]

  #transformation of data using rlog (which is suit for small datasets(n<30))
  dds_rlog <- rlog(dds_filtered)

  #meanSdPlot
  meanSd_plot <- meanSdPlot(assays(dds_rlog)[[1]], ranks=FALSE)

  #pcaPlot
  dds_pca = prcomp(assays(dds_rlog)[[1]])
  plot(dds_pca, type="l")
  p <- plotPCA(dds_rlog, intgroup="condition", returnData=TRUE)
  percentVar <- round(100*attr(p, "percentVar"))
  pca_plot <- ggplot(p, aes(PC1, PC2, color=condition, shape=condition, label=rownames(p))) +
    geom_point(size=2) +
    geom_text(hjust=0, vjust=1.6) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance"))


  meanSd_out <- file.path(dir, paste(outprefix, "_meanSd.pdf", sep = ""))
  pca_out <- file.path(dir, paste(outprefix, "_pca.pdf", sep = ""))

  pdf(meanSd_out, width = 15, height =10)
  print(meanSd_plot)
  dev.off()

  pdf(pca_out, width = 15, height =10)
  print(pca_plot)
  dev.off()
}

