#@param tool str, tools used to identify DEGs, edger or deseq
#@param top str, number of top significant genes to show
#@param colors list, not significant, up-regulated and down-regulated degs
#@param line bool, show threshold dashed lines or not
#@param gname, show gene names from, 0 from gene id column, 1 from gene annotation
#@param gnsep, the separator used in gene id column
#@param gncol, the column used to show in plot
#@param limit, the x limit for log2 fold change

library(ggplot2)
library(ggrepel)

get_degs_from_deseq <- function(treatment, control) {
	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC
	compare <- RNASUITE_COMPARE

	contrast <- c(compare, treatment, control)
	degs <- deseq_extract_degs(contrast)
	x <- as.data.frame(na.omit(degs))
	colnames(x)[which(names(x) == 'padj')] <- 'FDR'
	x$deg <- 'NS'
	x$deg[x$log2FoldChange >= logfc & x$FDR < fdr] <- 'Up'
	x$deg[x$log2FoldChange <= -logfc & x$FDR < fdr] <- 'Down'
	up_label <- paste('Up:', sum(x$deg == 'Up'))
	down_label <- paste('Down:', sum(x$deg == 'Down'))
	x$deg <- factor(x$deg, levels=c('Up', 'Down', 'NS'), labels=c(up_label, down_label, 'NS'))
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

rnasuite_degs_volcano_plot_update <- function(id=NULL, data=NULL, name=NULL, top=10,
	point_color=c('#ff006e', '#3a86ff', '#e9ecef'), show_vline=TRUE, vline_type='dashed', vline_color='black',
	vline_width=0.5, show_zline=FALSE, zline_type='dashed', zline_color='black', zline_width=0.5,
	show_hline=TRUE, hline_type='dashed', hline_color='black', hline_width=0.5, theme_name='bw',
	base_size=11, legend_position='top', x_limit=c(0, 0), point_size=0.5, ...) {

	fdr <- RNASUITE_FDR
	logfc <- RNASUITE_LOGFC

	if (!is.null(id)) {
		data <- rnasuite_get_data(id)
		name <- rnasuite_get_name(id)
	}

	if (is.null(data)) {
		return(NULL)
	}

	if (top > 0) {
		data$label <- rep('', nrow(data))
		data$label[order(data$FDR)[1:top]] <- data$gene[order(data$FDR)[1:top]]
	}

	p <- ggplot(data, aes(log2FoldChange, -log10(FDR), color=deg)) +
		geom_point(aes(fill=deg), size=point_size) +
		scale_fill_manual(values=point_color) +
		scale_color_manual(values=point_color)

	if (show_vline) {
		p <- p + geom_vline(xintercept=c(-logfc, logfc), linetype=vline_type, colour=vline_color, linewidth=vline_width)
	}

	if (show_zline) {
		p <- p + geom_vline(xintercept=0, linetype=zline_type, colour=zline_color, linewidth=zline_width)
	}

	if (show_hline) {
		p <- p + geom_hline(yintercept=-log10(fdr), linetype=hline_type, colour=hline_color, linewidth=hline_width)
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
		p <- p + geom_text_repel(aes(label=label), max.overlaps=Inf, key_glyph=draw_key_point)
	}

	if (!all(x_limit == 0)) {
		p <- p + lims(x=x_limit)
	}

	show(p)
	new <- as.integer(hgd_id()$id)
	rnasuite_put_plot(id, new, name, p, data)
	out <- list(c(1, name, new, 'deg_volcanoplot'))
	return(out)
}

rnasuite_degs_volcano_plot_run <- function(treatment, control, gname, gnsep, gncol, ...) {
	tool <- RNASUITE_DEGTOOL

	if (tool == 'deseq') {
		data <- get_degs_from_deseq(treatment, control)
	} else {
		data <- get_degs_from_edger(treatment, control)
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
	rnasuite_degs_volcano_plot_update(data=data, name=name)
}
