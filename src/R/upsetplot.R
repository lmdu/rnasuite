#@upsetplot_tool str, tools used to identify DEGs, edger or deseq
#@upsetplot_contrasts list, comparison contrasts
#@upsetplot_degtype, plot for degs, 0 for all, 1 for up-regulated, 2 for down-regulated

library(UpSetR)

get_degs_from_deseq <- function() {
	degs <- list()

	for (contrast in upsetplot_contrasts) {
		deseq_contrast[2:3] <<- contrast
		deseq_extract_degs()
		sig_degs <- deseq_sig_degs()
		label <- paste(contrast, collapse=' vs ')

		if (upsetplot_degtype == 1) {
			degs[[ label ]] = rownames(sig_degs[sig_degs$log2FoldChange>=deseq_logfc & sig_degs$padj<deseq_fdr, ])
		} else if (upsetplot_degtype == 2) {
			degs[[ label ]] = rownames(sig_degs[sig_degs$log2FoldChange<=-deseq_logfc & sig_degs$padj<deseq_fdr, ])
		} else {
			degs[[ label ]] = rownames(sig_degs[abs(sig_degs$log2FoldChange)>=deseq_logfc & sig_degs$padj<deseq_fdr, ])
		}
	}

	return(degs)
}

degs_upset_plot <- function() {
	if (upsetplot_tool == 'deseq') {
		data <- get_degs_from_deseq()
	}

	p <- upset(fromList(data),
		nsets = length(data),
		nintersects = 100,
		order.by = 'freq',
		point.size = 2.5,
		line.size = 0.8,
		text.scale = 1.8
	)

	show(p)
}
