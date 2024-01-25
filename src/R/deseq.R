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

	suppressWarnings(dds <- DESeqDataSetFromMatrix(
		countData = read_counts,
		colData = sample_info,
		design = as.formula(design_formula)
	))

	RNASUITE_DESEQ_RESULTS <<- DESeq(dds)
}

deseq_extract_degs <- function(contrast) {
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	degs <- results(RNASUITE_DESEQ_RESULTS,
		alpha = fdr,
		lfcThreshold = logfc,
		contrast = contrast
	)
	return(degs)
}

deseq_sig_degs <- function(degs) {
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	sig_degs <- na.omit(degs)
	sig_degs <- sig_degs[(sig_degs$padj < fdr & abs(sig_degs$log2FoldChange) >= logfc), ]
	sig_degs <- sig_degs[order(sig_degs$padj), ]
	sig_degs <- as.data.frame(sig_degs)
	return(sig_degs)
}

deseq_normalized_counts <- function() {
	norm_counts <- counts(RNASUITE_DESEQ_RESULTS, normalized=TRUE)
	norm_counts <- r_to_py(as.data.frame(norm_counts))
	return(norm_counts)
}

rnasuite_deseq_find_degs <- function(counts, samples, fdr, logfc, design, compare, treatment, control, ...) {
	RNASUITE_FDR <<- fdr
	RNASUITE_LOGFC <<- logfc
	RNASUITE_COMPARE <<- compare
	RNASUITE_TREATMENT <<- treatment
	RNASUITE_CONTROL <<- control
	RNASUITE_DEGTOOL <<- 'deseq'

	contrast <- c(compare, treatment, control)
	deseq_identify_degs(counts, samples, design)
	degs <- deseq_extract_degs(contrast)
	degs_list = r_to_py(deseq_sig_degs(degs))
	normal_count = deseq_normalized_counts()
	plotMA(degs)
	plot <- recordPlot()
	new <- as.integer(hgd_id()$id)
	name <- paste(treatment, 'vs', control, 'MA plot')
	old <- rnasuite_get_id(name)
	rnasuite_put_plot(old, new, name, plot, degs)

	out <- list(
		c(1, name, new, 'deseq_maplot'),
		c(0, 'Gene normalized read counts', normal_count, 'table'),
		c(0, paste(treatment, 'vs', control, 'DEGs list'), degs_list, 'table')
	)
	return(out)
}

rnasuite_deseq_extract_degs <- function(treatment, control, ...) {
	RNASUITE_TREATMENT <<- treatment
	RNASUITE_CONTROL <<- control

	compare <- RNASUITE_COMPARE
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	contrast <- c(compare, treatment, control)

	degs <- deseq_extract_degs(contrast)
	degs_list = r_to_py(deseq_sig_degs(degs))
	plotMA(degs)
	plot <- recordPlot()
	new <- as.integer(hgd_id()$id)
	name <- paste(treatment, 'vs', control, 'MA plot')
	old <- rnasuite_get_id(name)

	rnasuite_put_plot(old, new, name, plot, degs)
	out <- list(
		c(0, paste(treatment, 'vs', control, 'DEGs list'), degs_list, 'table'),
		c(1, name, new, 'deseq_maplot')
	)
	return(out)
}

rnasuite_deseq_ma_plot_update <- function(id=NULL, ...) {
	name <- rnasuite_get_name(id)
	data <- rnasuite_get_data(id)

	if (is.null(data)) {
		return(NULL)
	}

	plotMA(data, ...)
	plot <- recordPlot()
	new = as.integer(hgd_id()$id)
	rnasuite_put_plot(id, new, name, plot, data)
	out <- list(c(1, name, new, 'deseq_maplot'))
	return(out)
}
