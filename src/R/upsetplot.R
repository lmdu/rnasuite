#@upsetplot_tool str, tools used to identify DEGs, edger or deseq
#@upsetplot_contrasts list, comparison contrasts
#@upsetplot_degtype, plot for degs, 0 for all, 1 for up-regulated, 2 for down-regulated

library(UpSetR)

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

rnasuite_degs_upset_plot_update <- function(id=NULL, data=NULL, point_size=2.5, line_width=0.8, text_scale=1.8,
	order_by='freq', order_desc=TRUE, intersect_num=100, main_color='#3b3b3b', matrix_color='#3b3b3b',
	main_ylabel='Intersection Size', matrix_alpha=0.5, set_color='#3b3b3b', set_xlabel='Set Size',
	show_number='yes', number_angles=0, show_empty=NULL) {
	if (!is.null(id)) {
		name <- rnasuite_get_name(id)
		data <- rnasuite_get_data(id)
	} else {
		name <- "DEGs upset plot"
	}

	if (is.null(data)) {
		return(NULL)
	}

	p <- upset(fromList(data),
		nsets = length(data),
		nintersects = intersect_num,
		order.by = order_by,
		decreasing = order_desc,
		point.size = point_size,
		line.size = line_width,
		text.scale = text_scale,
		main.bar.color = main_color,
		mainbar.y.label = main_ylabel,
		matrix.color = matrix_color,
		matrix.dot.alpha = matrix_alpha,
		sets.bar.color = set_color,
		sets.x.label = set_xlabel,
		show.numbers = show_number,
		number.angles = number_angles,
		empty.intersections = show_empty
	)

	show(p)
	new <- as.integer(hgd_id()$id)
	rnasuite_put_plot(id, new, name, p, data)
	out <- list(c(1, name, new, 'deg_upsetplot'))
	return(out)
}

rnasuite_degs_upset_plot_run <- function(degtype, contrasts, ...) {
	tool <- RNASUITE_DEGTOOL

	if (tool == 'deseq') {
		data <- get_degs_from_deseq(degtype, contrasts)
	}

	rnasuite_degs_upset_plot_update(data=data)	
}
