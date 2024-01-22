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

get_degs_from_deseq <- function(fdr, logfc, compare, treatment, control) {
	contrast <- c(compare, treatment, control)
	deseq_extract_degs(fdr, logfc, contrast)
	x <- as.data.frame(na.omit(deseq_degs))
	colnames(x)[which(names(x) == 'padj')] <- 'FDR'
	x$deg <- 'NS'
	x$deg[x$log2FoldChange >= logfc & x$FDR < fdr] <- 'Up'
	x$deg[x$log2FoldChange <= -logfc & x$FDR < fdr] <- 'Down'
	up_label <- paste('Up:', sum(x$deg == 'Up'))
	down_label <- paste('Down:', sum(x$deg == 'Down'))
	x$deg <- factor(x$deg, levels=c('NS', 'Up', 'Down'), labels=c('NS', up_label, down_label))
	return(x)
}

get_degs_from_edger <- function(fdr, logfc) {
	x <- topTags(edger_degs, n=Inf)$table
	colnames(x)[which(names(x) == 'logFC')] <- 'log2FoldChange'
	x$deg <- 'NS'
	x$deg[x$log2FoldChange >= logfc & x$FDR < fdr] <- 'Up'
	x$deg[x$log2FoldChange <= -logfc & x$FDR < fdr] <- 'Down'
	up_label <- paste('Up:', sum(x$deg == 'Up'))
	down_label <- paste('Down:', sum(x$deg == 'Down'))
	x$deg <- factor(x$deg, levels=c('NS', 'Up', 'Down'), labels=c('NS', up_label, down_label))
	return(x)
}

rnasuite_degs_volcano_plot_update <- function(id=NULL, data=NULL, name=NULL, fdr=1, logfc=0.05, top=10,
	fill_colors=c(), label_colors=c(), show_vline=TRUE, vline_type='dash', vline_color='gray',
	vline_width=1, show_zline=FALSE, vline_type='dash', vline_color='black', vline_width=1,
	show_hline=TRUE, hline_type='dash', hline_color='gray', hline_width=1, theme_name='bw',
	base_size=11, legend_position='top', y_limit=c(0, 0)) {

	if (!is.null(id)) {
		data <- rnasuite_get_data(id)
		name <- rnasuite_get_name(id)
	}

	if (top > 0) {
		data$label <- rep('', nrow(data))
		data$label[order(data$FDR)[1:top]] <- data$gene[order(data$FDR)[1:top]]
	}

	p <- ggplot(data, aes(log2FoldChange, -log10(FDR), color=deg)) +
		geom_point(aes(fill=deg)) +
		scale_fill_manual(values=fill_colors) +
		scale_color_manual(values=label_colors) +

	if (show_vline) {
		p <- p + geom_vline(xintercept=c(-logfc, logfc), linetype=vline_type, colour=vline_color, size=vline_width)
	}

	if (show_zline) {
		p <- p + geom_vline(xintercept=0, linetype=zline_type, colour=zline_color, size=zline_width)
	}

	if (show_hline) {
		p <- p + geom_hline(yintercept=-log10(fdr), linetype=hline_type, colour=hline_color, size=hline_width)
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
	p <- p + theme(legend.position=legend_position, legend.title=element_blank())

	if (top > 0) {
		p <- p + geom_text_repel(aes(label=label), max.overlaps=100, key_glyph=draw_key_point)
	}

	if (!all(y_limit == 0)) {
		p <- p + lims(y=y_limit)
	}

	show(p)
	plot_id <- as.integer(hgd_id()$id)
	rnasuite_put_plot(id, plot_id, p, data)
	out <- list(c(1, name, plot_id, 'deg_distplot'))
	return(out)

	
}

rnasuite_degs_volcano_plot_show <- function(id) {
	p <- rnasuite_get_plot(id)
	show(p)
}

rnasuite_degs_volcano_plot_run <- function(tool, fdr, logfc, compare, treatment, control, gname, gnsep, gncol) {
	if (tool == 'deseq') {
		data <<- get_degs_from_deseq(fdr, logfc)
	} else {
		data <<- get_degs_from_edger(fdr, logfc)
	}

	if (gname == 0) {
		#data$label <- rep('', nrow(data))
		#data$label[order(data$FDR)[1:volcanoplot_top]] <- sapply(
		#	strsplit(rownames(data[order(data$FDR)[1:volcanoplot_top],]), volcanoplot_gnsep, fixed=T),
		#	function(x) {ifelse(is.na(x[volcanoplot_gncol]), x[1], x[volcanoplot_gncol])}
		#)
		data$gene <- sapply(
			strsplit(rownames(data), gnsep, fixed=TRUE),
			function(x) {ifelse(is.na(x[gncol]), x[1], x[gncol])}
		)
	}

	name = paste(treatment, 'vs', control, 'volcano plot')
	rnasuite_degs_volcano_plot_update(data=data, name=name, fdr=fdr, logfc=logfc)
}
