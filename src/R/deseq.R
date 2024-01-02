#@param dataframe read_counts, raw read count expression data
#@param dataframe sample_info, sample meta and group information
#@param str deseq_design, design formula for deg identification
#@param float deseq_fdr, FDR value for deg identification
#@param float deseq_logfc, log2foldchange for deg identification
#@param list deseq_contrast, results for which groups comparison
library(DESeq2)

deseq_identify_degs <- function() {
	#convert to integer
	read_counts <- round(read_counts)

	#remove genes with expression of 0
	read_counts <- read_counts[rowSums(read_counts[]) > 0, ]
	read_counts <- read_counts[, rownames(sample_info)]

	dds <- DESeqDataSetFromMatrix(
		countData = read_counts,
		colData = sample_info,
		design = as.formula(deseq_design)
	)

	deseq_results <<- DESeq(dds)
}

deseq_extract_degs <- function() {
	deseq_degs <<- results(deseq_results,
		alpha = deseq_fdr,
		lfcThreshold = deseq_logfc,
		contrast = deseq_contrast
	)
}

deseq_plot_degs <- function() {
	plotMA(deseq_degs)
}

deseq_sig_degs <- function() {
	sig_degs <- na.omit(deseq_degs)
	sig_degs <- sig_degs[(sig_degs$padj < deseq_fdr & abs(sig_degs$log2FoldChange) > deseq_logfc), ]
	sig_degs <- sig_degs[order(sig_degs$padj), ]
	sig_degs <- r_to_py(as.data.frame(sig_degs))
	return(sig_degs)
}

deseq_normalized_counts <- function() {
	norm_counts <- counts(deseq_results, normalized=TRUE)
	norm_counts <- r_to_py(as.data.frame(norm_counts))
	return(norm_counts)
}

deseq_return_degs <- function() {
	out <- list(
		normal_count = deseq_normalized_counts(),
		degs_list = deseq_sig_degs(),
		degs_versus = deseq_contrast
	)
	return(out)
}

deseq_analysis_pipeline <- function() {
	deseq_identify_degs()
	deseq_extract_degs()
	deseq_plot_degs()
	deseq_sig_degs()
	deseq_return_degs()
}

deseq_show_degs <- function() {
	deseq_extract_degs()
	deseq_plot_degs()
	out <- list(
		degs_list = deseq_sig_degs(),
		degs_versus = deseq_contrast
	)
	return(out)
}
