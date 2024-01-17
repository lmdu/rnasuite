#@param dataframe read_counts, raw read count expression data
#@param dataframe sample_info, sample meta and group information
#@param str deseq_design, design formula for deg identification
#@param float deseq_fdr, FDR value for deg identification
#@param float deseq_logfc, log2foldchange for deg identification
#@param list deseq_contrast, results for which groups comparison
library(DESeq2)

deseq_identify_degs <- function(read_counts, sample_info, design_formula) {
	#convert to dataframe
	read_counts <- py_to_r(r_to_py(read_counts))
	sample_info <- py_to_r(r_to_py(sample_info))

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
		alpha = fdr,
		lfcThreshold = logfc,
		contrast = contrast
	)
}

deseq_plot_degs <- function() {
	plotMA(deseq_degs)
	return(as.integer(hgd_id()$id))
}

deseq_sig_degs <- function(fdr, logfc) {
	sig_degs <- na.omit(deseq_degs)
	sig_degs <- sig_degs[(sig_degs$padj < fdr & abs(sig_degs$log2FoldChange) >= logfc), ]
	sig_degs <- sig_degs[order(sig_degs$padj), ]
	sig_degs <- as.data.frame(sig_degs)
	return(sig_degs)
}

deseq_normalized_counts <- function() {
	norm_counts <- counts(deseq_results, normalized=TRUE)
	norm_counts <- r_to_py(as.data.frame(norm_counts))
	return(norm_counts)
}

deseq_return_degs <- function(fdr, logfc, contrast) {
	normal_count = deseq_normalized_counts()
	degs_list = r_to_py(deseq_sig_degs(fdr, logfc))
	degs_plot = deseq_plot_degs()
	out <- list(
		c(0, 'Gene normalized read counts', normal_count, 'table'),
		c(0, paste(contrast[2], 'vs', contrast[3], 'DEGs'), degs_list, 'table'),
		c(1, paste(contrast[2], 'vs', contrast[3], 'MAplot'), degs_plot, 'deseq_maplot')
	)
	return(out)
}

rnasuite_deseq_find_degs <- function(counts, samples, fdr, logfc, design, compare, treatment, control, ...) {
	contrast <- c(compare, treatment, control)
	deseq_identify_degs(counts, samples, design)
	deseq_extract_degs(fdr, logfc, contrast)
	deseq_return_degs(fdr, logfc, contrast)
}

rnasuite_deseq_extract_degs <- function(fdr, logfc, compare, treatment, control, ...) {
	contrast <- c(compare, treatment, control)
	deseq_extract_degs(fdr, logfc, contrast)
	degs_list = r_to_py(deseq_sig_degs(fdr, logfc))
	degs_plot = deseq_plot_degs()
	out <- list(
		c(0, paste(treatment, 'vs', control, 'DEGs'), degs_list, 'table'),
		c(1, paste(treatment, 'vs', control, 'MAplot'), degs_plot, 'deseq_maplot')
	)
	return(out)
}

rnasuite_deseq_ma_plot_update <- function(treatment, control, ...) {
	plotMA(deseq_degs, ...)
	degs_plot = as.integer(hgd_id()$id)
	out <- list(c(1, paste(treatment, 'vs', control, 'MAplot'), degs_plot, 'deseq_maplot'))
	return(out)
}
