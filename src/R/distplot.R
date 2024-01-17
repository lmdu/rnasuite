#@distplot_tool str, tools used to identify DEGs, edger or deseq
#@distplot_contrasts list, comparison contrasts
#@distplot_type int, 0 for stacked bar plot, 1 for dodged bar plot
#@distplot_label bool, show value labels or not
#@distplot_rotate int, X labels rotate angle
#@distplot_colors list, color codes for up- and down-regulated bars

library(ggplot2)

deseq_degs_count <- function(fdr, logfc, compare, contrasts) {
	counts <- data.frame()
	groups <- c()

	for (contrast in contrasts) {
		comparison <- c(compare, contrast)
		deseq_extract_degs(fdr, logfc, comparison)
		sig_degs <- deseq_sig_degs(fdr, logfc)
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

edger_degs_count <- function() {

}

rnasuite_degs_dist_plot_update <- function(plot_type=0, show_label=FALSE, bar_colors=c('#E74C3C', '#3498DB'), x_rotate=0, plot_title=NULL,
	x_label=NULL, y_label='Counts', theme_name='bw', legend_title='DEGs') {
	counts <- distplot_data

	if (plot_type == 0) {
		counts$value <- as.numeric(counts$value) * c(1, -1)
		p <- ggplot(counts, aes(fill=condition, x=contrast, y=value)) +
			geom_bar(position='stack', stat='identity') +
			scale_y_continuous(labels = abs)

		if (show_label) {
			p <- p + geom_text(aes(label=abs(value), hjust=0.5, vjust=ifelse(value < 0, 1.3, -0.3)))
		}
	} else {
		counts$value <- as.numeric(counts$value)
		p <- ggplot(counts, aes(fill=condition, x=contrast, y=value)) +
			geom_bar(position='dodge', stat='identity')

		if (show_label) {
			p <- p + geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.3)
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

	p <- p + theme_func()

	if (x_rotate) {
		p <- p + theme(axis.text.x=element_text(angle = x_rotate, hjust=1))
	}

	show(p)
	return(as.integer(hgd_id()$id))
}

rnasuite_degs_dist_plot_run <- function(tool, fdr, logfc, compare, contrasts) {
	if (tool == 'deseq') {
		distplot_data <<- deseq_degs_count(fdr, logfc, compare, contrasts)
	} else if (tool == 'edger') {
		distplot_data <<- edger_degs_count()
	}

	rnasuite_degs_dist_plot_update()
}
