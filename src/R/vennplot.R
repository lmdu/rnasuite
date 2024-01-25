#@vennplot_tool str, tools used to identify DEGs, edger or deseq
#@vennplot_contrasts list, comparison contrasts
#@vennplot_percent bool, show percentage value labels or not
#@vennplot_colors list, color codes for each category
#@vennplot_opacity float, fill color alpha
#@vennplot_degtype, plot for degs, 0 for all, 1 for up-regulated, 2 for down-regulated

library(ggvenn)

get_degs_from_deseq <- function(degtype, contrasts) {
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	compare <- RNASUITE_COMPARE
	degs <- list()

	for (contrast in contrasts) {
		comparison <- c(compare, contrast)
		all_degs <- deseq_extract_degs(comparison)
		sig_degs <- deseq_sig_degs(all_degs)
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

rnasuite_degs_venn_plot_update <- function(id=NULL, data=NULL, fill_color=c('#E74C3C', '#3498DB', '#2ECC71', '#F1C40F'), ...) {
	if (!is.null(id)) {
		name <- rnasuite_get_name(id)
		data <- rnasuite_get_data(id)
	} else {
		name <- "DEGs venn plot"
	}

	if (is.null(data)) {
		return(NULL)
	}

	p <- ggvenn(data, fill_color = fill_color, ...)
	show(p)
	new <- as.integer(hgd_id()$id)
	rnasuite_put_plot(id, new, name, p, data)
	out <- list(c(1, name, new, 'deg_vennplot'))
	return(out)
}

rnasuite_degs_venn_plot_run <- function(degtype, contrasts, ...) {
	tool <- RNASUITE_DEGTOOL

	if (tool == 'deseq') {
		data <- get_degs_from_deseq(degtype, contrasts)
	}

	rnasuite_degs_venn_plot_update(data=data)
}
