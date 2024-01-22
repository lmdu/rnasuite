#@upsetplot_tool str, tools used to identify DEGs, edger or deseq
#@upsetplot_contrasts list, comparison contrasts
#@upsetplot_degtype, plot for degs, 0 for all, 1 for up-regulated, 2 for down-regulated

library(UpSetR)

get_degs_from_deseq <- function(fdr, logfc, degtype, compare, contrasts) {
	degs <- list()

	for (contrast in contrasts) {
		comparison <- c(compare, contrast)
		deseq_extract_degs(fdr, logfc, comparison)
		sig_degs <- deseq_sig_degs(fdr, logfc)
		label <- paste(contrast, collapse=' vs ')

		if (degtype == 1) {
			degs[[ label ]] = rownames(sig_degs[sig_degs$log2FoldChange>=logfc & sig_degs$padj<fdr, ])
		} else if (degtype == 2) {
			degs[[ label ]] = rownames(sig_degs[sig_degs$log2FoldChange<=-logfc & sig_degs$padj<fdr, ])
		} else {
			degs[[ label ]] = rownames(sig_degs[abs(sig_degs$log2FoldChange)>=logfc & sig_degs$padj<fdr, ])
		}
	}

	return(degs)
}

rnasuite_degs_upset_plot_update <- function() {
	p <- upset(fromList(upsetplot_data),
		nsets = length(upsetplot_data),
		nintersects = 100,
		order.by = 'freq',
		point.size = 2.5,
		line.size = 0.8,
		text.scale = 1.8
	)

	show(p)
}

rnasuite_degs_upset_plot_run <- function(tool, ...) {
	if (tool == 'deseq') {
		upsetplot_data <- get_degs_from_deseq(...)
	}

	rnasuite_degs_upset_plot_update()	
}
