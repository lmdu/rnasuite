#@vennplot_tool str, tools used to identify DEGs, edger or deseq
#@vennplot_contrasts list, comparison contrasts
#@vennplot_percent bool, show percentage value labels or not
#@vennplot_colors list, color codes for each category
#@vennplot_opacity float, fill color alpha
#@vennplot_degtype, plot for degs, 0 for all, 1 for up-regulated, 2 for down-regulated

library(ggvenn)

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

rnasuite_degs_venn_plot_update <- function(fill_color=c('#E74C3C', '#3498DB', '#2ECC71', '#F1C40F'), ...) {
	p <- ggvenn(vennplot_data, fill_color = fill_color, ...)
	show(p)
	plot_id <- as.integer(hgd_id()$id)
	out <- list(c(1, "DEGs venn plot", plot_id, 'deg_vennplot'))
	return(out)
}

rnasuite_degs_venn_plot_run <- function(tool, fdr, logfc, degtype, compare, contrasts) {
	if (tool == 'deseq') {
		vennplot_data <<- get_degs_from_deseq(fdr, logfc, degtype, compare, contrasts)
	}

	rnasuite_degs_venn_plot_update()
}
