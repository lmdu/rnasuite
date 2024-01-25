#@distplot_tool str, tools used to identify DEGs, edger or deseq
#@distplot_contrasts list, comparison contrasts
#@distplot_type int, 0 for stacked bar plot, 1 for dodged bar plot
#@distplot_label bool, show value labels or not
#@distplot_rotate int, X labels rotate angle
#@distplot_colors list, color codes for up- and down-regulated bars

library(ggplot2)

get_deseq_degs_count <- function(contrasts) {
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	compare <- RNASUITE_COMPARE
	counts <- data.frame()
	groups <- c()

	for (contrast in contrasts) {
		comparison <- c(compare, contrast)
		degs <- deseq_extract_degs(comparison)
		sig_degs <- deseq_sig_degs(degs)
		up_num <- sum(sig_degs$log2FoldChange>0)
		down_num <- sum(sig_degs$log2FoldChange<0)
		label <- paste(contrast, collapse=' vs ')
		counts <- rbind(counts, c(label, 'Up-regulated', up_num))
		counts <- rbind(counts, c(label, 'Down-regulated', down_num))
		groups <- append(groups, label)
	}

	colnames(counts) <- c('contrast', 'condition', 'value')
	counts$contrast <- factor(counts$contrast, levels=groups)
	counts$condition <- factor(counts$condition, levels=c('Up-regulated', 'Down-regulated'))
	return(counts)
}

get_edger_degs_count <- function(contrasts) {

}

rnasuite_degs_dist_plot_update <- function(id=NULL, data=NULL, plot_type=0, show_label=FALSE,
	bar_colors=c('#E74C3C', '#3498DB'), x_rotate=0, plot_title=NULL, x_label=NULL, y_label='Counts',
	theme_name='bw', legend_title='DEGs', base_size=11, ...) {

	if (!is.null(id)) {
		name <- rnasuite_get_name(id)
		data <- rnasuite_get_data(id)
	} else {
		name <- "DEGs distribution plot"
	}

	if (is.null(data)) {
		return(NULL)
	}

	counts <- data

	if (plot_type == 0) {
		counts$value <- as.numeric(counts$value) * c(1, -1)
		p <- ggplot(counts, aes(fill=condition, x=contrast, y=value)) +
			geom_bar(position='stack', stat='identity') +
			scale_y_continuous(labels = abs)

		if (show_label) {
			p <- p + geom_text(aes(label=abs(value), hjust=0.5, vjust=ifelse(value < 0, 1.3, -0.3)), size=base_size/3)
		}
	} else {
		counts$value <- as.numeric(counts$value)
		p <- ggplot(counts, aes(fill=condition, x=contrast, y=value)) +
			geom_bar(position='dodge', stat='identity')

		if (show_label) {
			p <- p + geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.3, size=base_size/3)
		}
	}

	p <- p + scale_fill_manual(values=bar_colors) +
		guides(fill=guide_legend(title=legend_title))

	if (length(plot_title)) {
		p <- p + ggtitle(plot_title)
	}

	if (length(x_label)) {
		p <- p + xlab(x_label)
	}

	if (length(y_label)) {
		p <- p + ylab(y_label)
	}

	theme_func <- switch(theme_name,
		'bw' = theme_bw,
		'classic' = theme_classic,
		'linedraw' = theme_linedraw,
		'minimal' = theme_minimal,
		'void' = theme_void,
		'light' = theme_light,
		'grey' = theme_grey,
		'gray' = theme_gray,
		'dark' = theme_dark
	)

	p <- p + theme_func(base_size=base_size)

	if (x_rotate) {
		p <- p + theme(axis.text.x=element_text(angle = x_rotate, hjust=1))
	}

	show(p)
	new <- as.integer(hgd_id()$id)
	rnasuite_put_plot(id, new, name, plot, data)
	out <- list(c(1, name, new, 'deg_distplot'))
	return(out)
}

rnasuite_degs_dist_plot_run <- function(contrasts, ...) {
	tool <- RNASUITE_DEGTOOL

	if (tool == 'deseq') {
		data <- get_deseq_degs_count(contrasts)
	} else if (tool == 'edger') {
		data <- get_edger_degs_count(contrasts)
	}

	rnasuite_degs_dist_plot_update(data=data)
}
