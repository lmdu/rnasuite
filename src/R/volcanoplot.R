#@volcanoplot_tool str, tools used to identify DEGs, edger or deseq
#@volcanoplot_top str, number of top significant genes to show
#@volcanoplot_colors list, not significant, up-regulated and down-regulated degs
#@volcanoplot_line bool, show threshold dashed lines or not
#@volcanoplot_gname, show gene names from, 0 from gene id column, 1 from gene annotation
#@volcanoplot_gnsep, the separator used in gene id column
#@volcanoplot_gncol, the column used to show in plot
#@volcanoplot_limit, the x limit for log2 fold change

library(ggplot2)
library(ggrepel)

get_degs_from_deseq <- function() {
	x <- as.data.frame(na.omit(deseq_degs))
	colnames(x)[which(names(x) == 'padj')] <- 'FDR'
	x$deg <- 'NS'
	x$deg[x$log2FoldChange >= deseq_logfc & x$FDR < deseq_fdr] <- 'Up'
	x$deg[x$log2FoldChange <= -deseq_logfc & x$FDR < deseq_fdr] <- 'Down'
	up_label <- paste('Up:', sum(x$deg == 'Up'))
	down_label <- paste('Down:', sum(x$deg == 'Down'))
	x$deg <- factor(x$deg, levels=c('NS', 'Up', 'Down'), labels=c('NS', up_label, down_label))
	return(x)
}

get_degs_from_edger <- function() {
	x <- topTags(edger_degs, n=Inf)$table
	colnames(x)[which(names(x) == 'logFC')] <- 'log2FoldChange'
	x$deg <- 'NS'
	x$deg[x$log2FoldChange >= edger_logfc & x$FDR < edger_fdr] <- 'Up'
	x$deg[x$log2FoldChange <= -edger_logfc & x$FDR < edger_fdr] <- 'Down'
	up_label <- paste('Up:', sum(x$deg == 'Up'))
	down_label <- paste('Down:', sum(x$deg == 'Down'))
	x$deg <- factor(x$deg, levels=c('NS', 'Up', 'Down'), labels=c('NS', up_label, down_label))
	return(x)
}

degs_volcano_plot <- function() {
	if (volcanoplot_tool == 'deseq') {
		data = get_degs_from_deseq()
		logfc <- deseq_logfc
		fdr <- deseq_fdr
	} else {
		data = get_degs_from_edger()
		logfc <- edger_logfc
		fdr <- edger_fdr
	}

	if (volcanoplot_gname == 0) {
		data$label <- rep('', nrow(data))
		data$label[order(data$FDR)[1:volcanoplot_top]] <- sapply(
			strsplit(rownames(data[order(data$FDR)[1:volcanoplot_top],]), volcanoplot_gnsep, fixed=T),
			function(x) {ifelse(is.na(x[volcanoplot_gncol]), x[1], x[volcanoplot_gncol])}
		)
	}

	p <- ggplot(data, aes(log2FoldChange, -log10(FDR), color=deg)) +
		geom_point(aes(fill=deg)) +
		scale_fill_manual(values=volcanoplot_colors) +
		scale_color_manual(values=volcanoplot_colors) +
		geom_vline(xintercept=c(-logfc, logfc), linetype='dashed') +
		geom_hline(yintercept=-log10(fdr), linetype='dashed') +
		theme_bw() +
		theme(legend.position="top", legend.title=element_blank())

	if (volcanoplot_top > 0) {
		p <- p + geom_text_repel(aes(label=label), max.overlaps=100, key_glyph=draw_key_point)
	}

	if (volcanoplot_limit > 0) {
		p <- p + xlim(-volcanoplot_limit, volcanoplot_limit)
	}

	show(p)
}
