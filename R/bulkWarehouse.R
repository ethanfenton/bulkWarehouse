#' Multiple Testing Recorrection of DESeq2 results
#'
#' After identifying genes that should be filtered out of a bulk analysis (e.g. due to low read counts), recorrect the pvalues with multiple testing correction via the Benjamini Hochberg procedure
#' @param res.obj a results dataframe from DEseq2
#' @param filt.out.genes.vec A vector of gene ensembl ids to remove from the anlaysis for the purposed of multiple test correction
#' @return A dataframe same as the input dataframe but with two extra columns. pvalue.filt now contains the original pvalue or NA if the gene was filtered out. padj_BH now contains the BenjaminiHochberg adjusted pvalue based on the genes that pass the filter
#' @examples
#' recorrected_DESeq2_resutls <- readjust_DESeq2_results(res.obj=my_DESeq2_results, filt.out.genes.vec=c("ENSMUSG00000000766","ENSMUST00000036653.5"));
#' @export
readjust_DESeq2_results <- function(res.obj, filt.out.genes.vec){
  res.obj$pvalue.filt <- res.obj$pvalue
  res.obj$pvalue.filt[which(res.obj$ensembl_id %in% filt.out.genes.vec)] <- NA
  #res.obj$pvalue.filt[ filt.out.genes.vec ] <- NA
  res.obj$padj_BH <- p.adjust(res.obj$pvalue.filt, method="BH")
  return(res.obj)
}

#' Permutation test for Bulk RNAseq
#'
#' Generate an empirical p-value to determine if a given gene set is preferentially represented, more frequently than chance, in the top of a differnetilaly expressed gene list
#' @param DE_list A df of differentially expressed genes from DEseq2 with a pvalue and/or padj and/or log2FoldChange, IMPORTANT: the name of the gene must be in a column titled Gene.name and be in the proper format e.g. Tbr2
#' @param test_genes A list of gene names that you wish to test for overrepresentation at the top of your DE gene df
#' @param pval A pvalue cutoff for which genes with a smaller pvalue will be allowed to be considered part of the "top DE genes"
#' @param l2fc A log2FoldChange cutoff for which genes with a grater absolute value of log2FoldChange will be allowed to be considered part of the "top DE genes"
#' @param adjusted A boolean that determines if the pvalue or adjusted pvalue will be considered with the above pval parameter
#' @param runs The number of simulations to run
#' @param binwidth Binwodth parameter for plotting resulting emperical pvalue histogram
#' @return Empirical pvalue
#' @examples
#' pvalue <- perm_test(DE_list=DEseq2_results, test_genes=toupper(my_gwas_genes$Gene.name), pval=0.05, runs=10000, l2fc=0.3);
#' @export
perm_test <- function(DE_list, test_genes, pval, l2fc=0, adjusted=FALSE, runs=10000, binwidth=0.1){
  if(adjusted){
    DE_filtered <- DE_list %>% filter(padj <= pval) %>% filter(abs(log2FoldChange) > l2fc)
  }else{
    DE_filtered <- DE_list %>% filter(pvalue <= pval) %>% filter(abs(log2FoldChange) > l2fc)
  }
  test_ratio <- sum(toupper(test_genes) %in% toupper(DE_filtered$Gene.name))/length(DE_filtered$Gene.name)
  print(paste("Test Ratio: ", test_ratio))

  random_ratios <- c()
  for(i in 1:runs){
    random_genes <- sample(DE_list$Gene.name, size=length(test_genes), replace=FALSE)
    #print(random_genes[27])
    random_ratio <- sum(toupper(random_genes) %in% toupper(DE_filtered$Gene.name))/length(DE_filtered$Gene.name)
    #print(random_ratio)
    random_ratios[length(random_ratios)+1] <- random_ratio
  }
  random_ratios_df <- data.frame(the_random_ratios=random_ratios*100)
  p <- ggplot(random_ratios_df) +
    geom_histogram(aes(x=the_random_ratios), binwidth=0.1) +
    geom_vline(aes(xintercept=gwas_ratio*100), color="red") +
    geom_vline(aes(xintercept=quantile(random_ratios, 0.95)*100), color="green") +
    geom_vline(aes(xintercept=quantile(random_ratios, 0.05)*100), color="green") +
    ggtitle("Permutation Distribution,\n red is our statistic and\n green is p=0.05 or p=0.95")
  plot(p)
  print(
    paste( "Your pvalue is:",( 1-ecdf(random_ratios)(test_ratio) ) )
  )
  return(1-ecdf(random_ratios)(test_ratio))
}
