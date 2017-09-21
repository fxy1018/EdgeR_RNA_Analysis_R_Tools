setwd("~/RNA_Seq/RNA_Seq/IWP0004JJ/Quantifier")

dir <- "~/RNA_Seq/RNA_Seq/IWP0004JJ/Quantifier"
files <- c("S1_sequence_non_rRNA_count.txt","S2_sequence_non_rRNA_count.txt","S3_sequence_non_rRNA_count.txt",
           "S4_sequence_trimmed_non_rRNA_count.txt","S5_sequence_trimmed_non_rRNA_count.txt","S6_sequence_trimmed_non_rRNA_count.txt",
           "S7_sequence_non_rRNA_count.txt","S8_sequence_non_rRNA_count.txt")

files <- paste(files, ".deseq2", sep="")
condition <- factor(c(0,0,1,1,2,2,3,3))
outprefix <- "IWP0004JJ"
qcPlots(files, dir, condition, outprefix)


files <- c("S1_sequence_non_rRNA_count.txt","S2_sequence_non_rRNA_count.txt","S3_sequence_non_rRNA_count.txt",
           "S4_sequence_trimmed_non_rRNA_count.txt","S5_sequence_trimmed_non_rRNA_count.txt","S6_sequence_trimmed_non_rRNA_count.txt",
           "S7_sequence_non_rRNA_count.txt","S8_sequence_non_rRNA_count.txt")
files <- paste(files, ".deseq2", sep="")
dir <- "~/RNA_Seq/RNA_Seq/IWP0004JJ/Quantifier"
group <- factor(c(0,0,1,1,2,2,3,3))

fit <- edgeRAnaFit(files, dir, group, outprefix)


coefs = list(2,3,4,c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))
comparisons = list("condition_1_vs_0","condition_2_vs_0","condition_3_vs_0",
             "condition_2_vs_1","condition_3_vs_1","condition_3_vs_2")
dfGenes <- list()
lrts <- list()
for (i in (1:length(coefs))){
  coef = coefs[[i]]
  return = edgeRdfGenes(fit, coef)
  lrts[[i]] = return[[1]]
  dfGenes[[i]] = return[[2]]
}

df = edgeRAnaWrite(dfGenes,comparisons, dir, outprefix)

#pathwway analysis
#kegg pathway
pathways = edgeRAnaKEGG(lrts, dfCutOff=0.05, dir, comparisons, outprefix)

#get significant pathways
sign_path <- getPathway(pathways, fdr=0.05, out= T, dir, outprefix)

#map genes to significant pathway
mapped_genes <- edgeRAnaMapedGene2Pathway(sign_path,df, out=T, dir ,outprefix)

#GO
gos= edgeRAnaGO(lrts, fdr=0.05, dir, comparisons, outprefix)


#reactome pathway
reactome <- edgeRAnaReactome(df, fdr=0.05, dir, comparisons, outprefix)

#pathview datat
edgeRAnaPathview(mapped_genes, "~/RNA_Seq/RNA_Seq/IWP0004JJ/Quantifier/pathview")


