#@vennplot_tool str, tools used to identify DEGs, edger or deseq
#@vennplot_contrasts list, comparison contrasts
#@vennplot_percent bool, show percentage value labels or not
#@vennplot_colors list, color codes for each category
#@vennplot_opacity float, fill color alpha
#@vennplot_degtype, plot for degs, 0 for all, 1 for up-regulated, 2 for down-regulated

library(ggvenn)

get_degs_from_deseq <- function() {
	degs <- list()

	for (contrast in vennplot_contrasts) {
		deseq_contrast[2:3] <<- contrast
		deseq_extract_degs()
		sig_degs <- deseq_sig_degs()
		label <- paste(contrast, collapse=' vs ')

		if (vennplot_degtype == 1) {
			degs[[ label ]] = rownames(sig_degs[sig_degs$log2FoldChange>=deseq_logfc & sig_degs$padj<deseq_fdr, ])
		} else if (vennplot_degtype == 2) {
			degs[[ label ]] = rownames(sig_degs[sig_degs$log2FoldChange<=-deseq_logfc & sig_degs$padj<deseq_fdr, ])
		} else {
			degs[[ label ]] = rownames(sig_degs[abs(sig_degs$log2FoldChange)>=deseq_logfc & sig_degs$padj<deseq_fdr, ])
		}
	}

	return(degs)
}

degs_venn_plot <- function() {
	if (vennplot_tool == 'deseq') {
		data <- get_degs_from_deseq()
	}

	p <- ggvenn(data, digits = 2,
		fill_color = vennplot_colors,
		fill_alpha = vennplot_opacity,
		show_percentage = vennplot_percent
	)

	show(p)
}
