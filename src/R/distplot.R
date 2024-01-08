#@distplot_tool str, tools used to identify DEGs, edger or deseq
#@distplot_contrasts list, comparison contrasts
#@distplot_type int, 0 for stacked bar plot, 1 for dodged bar plot
#@distplot_label bool, show value labels or not
#@distplot_rotate int, X labels rotate angle
#@distplot_colors list, color codes for up- and down-regulated bars

library(ggplot2)

deseq_degs_count <- function() {
	counts <- data.frame()
	groups <- c()

	for (contrast in distplot_contrasts) {
		deseq_contrast[2:3] <<- contrast
		deseq_extract_degs()
		sig_degs <- deseq_sig_degs()
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

degs_dist_plot <- function() {
	if (distplot_tool == 'deseq') {
		counts <- deseq_degs_count()
	} else if (distplot_tool == 'edger') {
		counts <- edger_degs_count()
	}

	if (distplot_type == 0) {
		counts$value <- as.numeric(counts$value) * c(1, -1)
		p <- ggplot(counts, aes(fill=condition, x=contrast, y=value)) +
			geom_bar(position='stack', stat='identity') +
			scale_y_continuous(labels = abs)

		if (distplot_label) {
			p <- p + geom_text(aes(label=abs(value), hjust=0.5, vjust=ifelse(value < 0, 1.3, -0.3)))
		}
	} else {
		counts$value <- as.numeric(counts$value)
		p <- ggplot(counts, aes(fill=condition, x=contrast, y=value)) +
			geom_bar(position='dodge', stat='identity')

		if (distplot_label) {
			p <- p + geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.3)
		}
	}

	p <- p + scale_fill_manual(values=distplot_colors) +
		guides(fill=guide_legend(title="DEGs")) +
		theme_bw() +
		xlab('') +
		ylab('Counts')

	if (distplot_rotate) {
		p <- p + theme(axis.text.x=element_text(angle = distplot_rotate, hjust=1))
	}

	show(p)
	return(NULL)
}
