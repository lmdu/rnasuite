#@param dataframe read_counts, raw read count expression data
#@param dataframe sample_info, sample meta and group information
#@param str deseq_design, design formula for deg identification
#@param float deseq_fdr, FDR value for deg identification
#@param float deseq_logfc, log2foldchange for deg identification
#@param list deseq_contrast, results for which groups comparison
library(DESeq2)

DeseqIdentifyDegs <- function(read.counts, sample.info, design.formula) {
	#convert to dataframe
	read.counts <- RnasuitePandasToDataframe(read.counts)
	sample.info <- RnasuitePandasToDataframe(sample.info)

	#convert to integer
	read.counts <- round(read.counts)

	#remove genes with expression of 0
	read.counts <- read.counts[rowSums(read.counts[]) > 0, ]
	read.counts <- read.counts[, rownames(sample.info)]

	suppressWarnings(dds <- DESeqDataSetFromMatrix(
		countData = read.counts,
		colData = sample.info,
		design = as.formula(design.formula)
	))

	RNASUITE.DESEQ.RESULTS <<- DESeq(dds)
}

DeseqExtractDegs <- function(fdr, logfc, contrast) {
	degs <- results(RNASUITE.DESEQ.RESULTS,
		alpha = fdr,
		lfcThreshold = logfc,
		contrast = contrast
	)
	return(degs)
}

DeseqGetSigificantDegs <- function(fdr, logfc, degs) {
	sig.degs <- na.omit(degs)
	sig.degs <- sig.degs[(sig.degs$padj < fdr & abs(sig.degs$log2FoldChange) >= logfc), ]
	sig.degs <- sig.degs[order(sig.degs$padj), ]
	sig.degs <- as.data.frame(sig.degs)
	return(sig.degs)
}

DeseqGetNormalizedCounts <- function() {
	norm.counts <- counts(RNASUITE.DESEQ.RESULTS, normalized=TRUE)
	norm.counts <- as.data.frame(norm.counts)
	norm.counts <- RnasuiteDataframeToPandas(norm.counts)
	return(norm.counts)
}

Rnasuite_deseq_extract_degs <- function(treatment, control, ...) {
	RNASUITE_TREATMENT <<- treatment
	RNASUITE_CONTROL <<- control

	compare <- RNASUITE_COMPARE
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	contrast <- c(compare, treatment, control)

	degs <- DeseqExtractDegs(fdr, logfc, contrast)
	degs_list = r_to_py(deseq_sig_degs(degs))
	plotMA(degs)
	plot <- recordPlot()
	new <- as.integer(unigd::ugd_id()$id)
	name <- paste(treatment, 'vs', control, 'MA plot')
	old <- rnasuite_get_id(name)

	rnasuite_put_plot(old, new, name, plot, degs)
	out <- list(
		c(0, paste(treatment, 'vs', control, 'DEGs list'), degs_list, 'table'),
		c(1, name, new, 'deseq_maplot')
	)
	return(out)
}

RnasuiteDeseqFindDegs <- function(counts=NULL, samples=NULL, design=NULL, fdr,
	logfc, compare, treatment, control, ...) {
	contrast <- c(compare, treatment, control)
	result <- list()

	if (!is.null(counts)) {
		DeseqIdentifyDegs(counts, samples, design)
		normal.counts <- DeseqGetNormalizedCounts()
		item <- list(type=0, name='Normalized read counts', data=normal.counts)
		result <- append(result, item)
	}

	degs <- DeseqExtractDegs(fdr, logfc, contrast)
	degs.list <- DeseqGetSigificantDegs(fdr, logfc, degs)
	degs.list = RnasuiteDataframeToPandas(degs.list)

	plotMA(degs)
	plot <- recordPlot()
	plot.id <- RnasuiteGetPlotCurrentId()
	plot.name <- paste(treatment, 'vs', control, 'MA plot')
	old.id <- RnasuiteGetPlotIdByName(plot.name)

	RnasuiteSavePlot(old.id, plot.id, plot.name, degs)
	item <- list(type=1, name=plot.name, id=plot.id, plot='deseq_maplot')
	result <- append(result, item)
	item <- list(type=0, name=paste(treatment, 'vs', control, 'DEGs list'), data=degs.list)
	result <- append(result, item)

	return(result)
}

rnasuite_deseq_ma_plot_update <- function(id=NULL, ...) {
	name <- rnasuite_get_name(id)
	data <- rnasuite_get_data(id)

	if (is.null(data)) {
		return(NULL)
	}

	plotMA(data, ...)
	plot <- recordPlot()
	new = as.integer(unigd::ugd_id()$id)
	rnasuite_put_plot(id, new, name, plot, data)
	out <- list(c(1, name, new, 'deseq_maplot'))
	return(out)
}
