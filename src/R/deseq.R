#@param dataframe read_counts, raw read count expression data
#@param dataframe sample_info, sample meta and group information
#@param str deseq_design, design formula for deg identification
#@param float deseq_fdr, FDR value for deg identification
#@param float deseq_logfc, log2foldchange for deg identification
#@param list deseq_contrast, results for which groups comparison
library(DESeq2)

deseq_identify_degs <- function(read_counts, sample_info, design_formula) {
	#convert to integer
	read_counts <- round(read_counts)

	#remove genes with expression of 0
	read_counts <- read_counts[rowSums(read_counts[]) > 0, ]
	read_counts <- read_counts[, rownames(sample_info)]

	dds <- DESeqDataSetFromMatrix(
		countData = read_counts,
		colData = sample_info,
		design = as.formula(design_formula)
	)

	deseq_results <<- DESeq(dds)
}

deseq_extract_degs <- function(fdr, logfc, contrast) {
	deseq_degs <<- results(deseq_results,
		alpha = deseq_fdr,
		lfcThreshold = deseq_logfc,
		contrast = contrast
	)
}

deseq_plot_degs <- function() {
	plotMA(deseq_degs)
	return(as.integer(hgd_id()$id))
}

deseq_sig_degs <- function(fdr, logfc) {
	sig_degs <- na.omit(deseq_degs)
	sig_degs <- sig_degs[(sig_degs$padj < deseq_fdr & abs(sig_degs$log2FoldChange) >= deseq_logfc), ]
	sig_degs <- sig_degs[order(sig_degs$padj), ]
	sig_degs <- as.data.frame(sig_degs)
	return(sig_degs)
}

deseq_normalized_counts <- function() {
	norm_counts <- counts(deseq_results, normalized=TRUE)
	norm_counts <- r_to_py(as.data.frame(norm_counts))
	return(norm_counts)
}

deseq_return_degs <- function(contrast) {
	out <- list(
		normal_count = deseq_normalized_counts(),
		degs_list = r_to_py(deseq_sig_degs()),
		degs_plot = deseq_plot_degs(),
		plot_type = "deseq_maplot",
		degs_versus = contrast
	)
	return(out)
}

rnasuite_deseq_find_degs <- function(counts, samples, fdr, logfc, design, compare, treatment, control) {
	contrast <- c(compare, treatment, control)
	deseq_identify_degs(counts, samples, design)
	deseq_extract_degs(fdr, logfc, contrast)
	deseq_return_degs(contrast)
}

rnasuite_deseq_extract_degs <- function(fdr, logfc, compare, treatment, control) {
	contrast <- c(compare, treatment, control)
	deseq_extract_degs(fdr, logfc, contrast)
	out <- list(
		degs_list = r_to_py(deseq_sig_degs()),
		degs_plot = deseq_plot_degs(),
		plot_type = "deseq_maplot",
		degs_versus = contrast
	)
	return(out)
}

rnasuite_deseq_ma_plot <- function(...) {
	plotMA(deseq_degs, ...)
	return(as.integer(hgd_id()$id))
}
