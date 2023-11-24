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
#' @param de_df A df of differentially expressed genes from DEseq2 with a pvalue and/or padj and/or log2FoldChange, IMPORTANT: the name of the gene must be in a column titled Gene.name and be in the proper format e.g. Tbr2
#' @param test_df A df of genes that you wish to test for overrepresentation at the top of your DE gene df
#' @param de_pval A pvalue cutoff for which genes from the de_df with a smaller pvalue will be allowed to be considered part of the "top DE genes"
#' @param de_l2fc_min A log2FoldChange cutoff for which genes from the de_df with a greater absolute value of log2FoldChange will be allowed to be considered part of the "top DE genes"
#' @param test_pval The pval cutoff used to select genes for inclusion from the test_df in whatever units the test_df_pvar is in e.g. pvalue or fdr
#' @param test_fc_min The absolute-value-of-fold-change cutoff used to select genes for inclusion from the test_df in whatever units the test_fc_min is in e.g. log2 or non-log
#' @param take_top_test The number of genes of which to limit the test set
#' @param take_top_de The number of genes of which to limit the top DE set
#' @param adjusted boolean that specifies to use pvalue or padj for the top DE (assumes DESeq2 formatted df)
#' @param runs The number of simulations to run
#' @param my_binwidth Binwidth parameter for plotting resulting emperical enrichment histogram
#' @param verts The quantile to plot the vertical bars to help vizualize significance. If this is 5 then the 5th and 95th quantiles are plotted, else the 0.25 and 0.975 are plotted
#' @param reverse A boolean when FALSE: (num top_test_genes in top_de_genes)/(num top_de_genes) and when TRUE: (num top_test_genes in top_de_genes)/(num top_test_genes)
#' @param test_df_pvar The variable within the test_df that stores the pvalue for filtering
#' @param test_df_fc_var The variable within the test_df that stores the fold change metric for filtering
#' @return A named list with the item called plot being a ggplot object, and the item called emp_pval being the emperical pvalue as determined by the quantile in the distribution that the enrichment value for the test_df gene set falls
#' @examples perm_test(de_df=res.sal.lsd.pfilt, test_df=twas_input, de_pval=0.1, de_l2fc_min=0, test_pval=0.05, test_fc_min=0, take_top_test=500, take_top_de=500, adjusted=FALSE, runs=100000, my_binwidth=0.2, verts="5", reverse=FALSE, test_df_pvar="SCZ.fdr", test_df_fc_var="SCZ.log2FC")
#' @export

perm_test <- function(de_df=res.sal.lsd.pfilt, test_df=test_scz_asd_bp_Geschwind_2018_filt_SCZ,
                       de_pval=0.1, de_l2fc_min=0,
                       test_pval=0.005, test_fc_min=0,
                       take_top_test=10000, take_top_de=10000,
                       adjusted=FALSE, runs=10000, my_binwidth=0.2, verts="fif", reverse=FALSE,
                       test_df_pvar="SCZ.fdr", test_df_fc_var="SCZ.log2FC"){

  # Capture the arguments to put as plot title for record keeping
  args_list <- as.list(match.call()[-1])
  args_strings <- sapply(args_list, deparse)
  title_string <- paste(names(args_strings), args_strings, sep = "=", collapse = " ; ")
  #print(title_string)


  top_test_df <- filter(test_df, get(test_df_pvar) < test_pval & get(test_df_fc_var) > test_fc_min & gene_name %in% toupper(de_df$Gene.name)) %>% arrange(desc(abs(get(test_df_fc_var)))) %>% head(take_top_test)
  top_test_genes <- top_test_df$gene_name
  print(paste("Number of test Genes used as input:" ,length(top_test_genes)))
  #print(paste("the VERTS:",verts))

  if(adjusted){
    top_de_df <- de_df %>% filter(padj.filt <= de_pval) %>% filter(abs(log2FoldChange) > de_l2fc_min) %>% arrange(desc(abs(log2FoldChange))) %>% head(take_top_de)
  }else{
    top_de_df <- de_df %>% filter(pvalue.filt <= de_pval) %>% filter(abs(log2FoldChange) > de_l2fc_min) %>% arrange(desc(abs(log2FoldChange))) %>% head(take_top_de)
  }
  print(paste("Bulk set size:", nrow(de_df)))
  print(paste("Top DE genes in our Bulk set size:", nrow(top_de_df)))
  if(reverse){
    test_ratio <- sum(toupper(top_de_df$Gene.name) %in% toupper(top_test_genes)) / length(top_test_genes)
  }else{
    test_ratio <- sum(toupper(top_test_genes) %in% toupper(top_de_df$Gene.name)) / length(top_de_df$Gene.name)
  }
  print(paste("User Geneset Ratio: ", test_ratio))
  print("Intersection Genes:")
  #print(top_test_genes[toupper(top_test_genes) %in% toupper(top_de_df$Gene.name)])))
  intersect_genes <- top_test_genes[which(top_test_genes %in% toupper(top_de_df$Gene.name))]
  intersect_genes_string <- paste(intersect_genes, collapse = ", ")
  title_string <- paste(title_string,"\nIntersection Genes:", intersect_genes_string)

  random_ratios <- c()
  for(i in 1:runs){
    cat("\r","Iterations:", i)
    if(reverse){
      #my genes in test
      random_genes <- sample(de_df$Gene.name, size=length(top_de_genes), replace=FALSE)
      random_ratio <- sum(toupper(random_genes) %in% toupper(top_de_df$Gene.name))/length(top_test_genes)
    }else{
      #test in my deg
      random_genes <- sample(de_df$Gene.name, size=length(top_test_genes), replace=FALSE)
      random_ratio <- sum(toupper(random_genes) %in% toupper(top_de_df$Gene.name))/length(top_de_df$Gene.name)
    }
    #print(random_ratio)
    random_ratios[length(random_ratios)+1] <- random_ratio
  }
  random_ratios_df <- data.frame(the_random_ratios=random_ratios*100)
  #random_ratios_df <- na.omit(random_ratios_df)
  p <- NULL
  if(verts=="fif"){
    p <- ggplot(random_ratios_df) +
      geom_histogram(aes(x=the_random_ratios), binwidth =my_binwidth) +
      geom_vline(aes(xintercept=test_ratio*100), color="orange",size=1.5) +
      geom_vline(aes(xintercept=quantile(random_ratios, 0.95)*100), color="blue", linetype="dashed",size=1.5) +
      geom_vline(aes(xintercept=quantile(random_ratios, 0.05)*100), color="blue", linetype="dashed",size=1.5) +
      xlab("Prevalence of Ratio") +
      ylab("Frequency of Result") +
      ggtitle(title_string) +
      theme_bw()
  }else{
    p <- ggplot(random_ratios_df) +
      geom_histogram(aes(x=the_random_ratios), binwidth =my_binwidth) +
      geom_vline(aes(xintercept=test_ratio*100), color="orange",size=1.5) +
      geom_vline(aes(xintercept=quantile(random_ratios, 0.975)*100), color="pink", linetype="dashed",size=1.5) +
      geom_vline(aes(xintercept=quantile(random_ratios, 0.025)*100), color="pink", linetype="dashed",size=1.5) +
      xlab("Prevalence Ratio") +
      ylab("Frequency of Result") +
      ggtitle(title_string) +
      theme_bw()
  }
  emp.pval = 1-(ecdf(random_ratios)(test_ratio))
  print(
    paste( "Your pvalue is:",( 1-(ecdf(random_ratios)(test_ratio)) ) )
  )
  p <- p + theme(
    axis.text.x = element_text(size = 50),
    axis.text.y = element_text(size = 50),
    axis.title.x = element_text(size = 50),
    axis.title.y = element_text(size = 50),
    plot.title = element_text(size = 5)
  )

  plot(p, main=title_string)

  return(list(plot = p, emp_pval = emp.pval))
}

#' ampSearch C++ Function
#'
#' This function finds matching sequences in two vectors of sequences. This was made for the purpose of finding primer and cell barcodes in amplicon sequencing data
#'
#' @param short_strings A CharacterVector of strings for which to search (e.g. cell barcodes)
#' @param long_strings A CharacterVector of strings in which to search (e.g. amplicon seqeunces)
#' @param show_progress A Bool determining if the percent progress should be output
#' @return A 2D CharacterMatrix with a long_string in column1 and a matching short_string in column2, with one entry per column per row. A long string can have multiple short string matches (theoreticall) and a short string can have multiple long string matches (more likely scenario)
#' @export
#' @useDynLib bulkWarehouse
#' @importFrom Rcpp sourceCpp
#' @examples
#' ampSearch(short_strings=cell_barcodes_20000x_16bp, long_strings=amplicon_sequences_50000x_350bp, show_progres=True)
ampSearch <- function(short_strins, long_strings, show_progress=False) {
  .Call('_bulkWarehouse_ampSearch', PACKAGE = 'bulkWarehouse',
        short_strings = short_strings, 
        long_strings = long_strings, 
        show_progress = show_progress)
}
