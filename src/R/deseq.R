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
		result[[length(result)+1]] <- list(
			type = 'table',
			name = 'Normalized read counts',
			data = normal.counts,
			kind = 'normal_count'
		)
	}

	degs <- DeseqExtractDegs(fdr, logfc, contrast)
	degs.list <- DeseqGetSigificantDegs(fdr, logfc, degs)
	degs.list = RnasuiteDataframeToPandas(degs.list)
	result[[length(result)+1]] <- list(
		type = 'table',
		name = paste(treatment, 'vs', control, 'DEGs list'),
		data = degs.list,
		kind = 'degs_list'
	)

	plotMA(degs)
	plot <- recordPlot()
	plot.id <- RnasuiteGetPlotCurrentId()
	plot.name <- paste(treatment, 'vs', control, 'MA plot')
	old.id <- RnasuiteGetPlotIdByName(plot.name)
	RnasuiteSavePlot(old.id, plot.id, plot.name, plot, degs)
	result[[length(result)+1]] <- list(
		type = 'plot',
		name = plot.name,
		data = plot.id,
		kind = 'deseq_maplot'
	)

	return(result)
}

RnasuiteDeseqMaplotUpdate <- function(id=NULL, ...) {
	chart = RnasuiteGetPlot(id)

	if (is.null(chart$data)) {
		return(NULL)
	}

	plotMA(chart$data, ...)
	plot <- recordPlot()
	plot.id <- RnasuiteGetPlotCurrentId()
	RnasuiteSavePlot(chart$id, plot.id, chart$name, plot, chart$data)

	out <- list(list(
		type = 'plot',
		name = chart$name,
		data = plot.id,
		kind = chart$code
	))
	return(out)
}
