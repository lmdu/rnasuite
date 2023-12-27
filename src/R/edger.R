#@param dataframe read_counts, raw read count expression data
#@param dataframe sample_info, sample meta and group information
#@param str edger_design, design formula for deg identification
#@param float edger_fdr, FDR value for deg identification
#@param float edger_logfc, log2foldchange for deg identification
#@param str edger_group, find DEGs between given group
#@param int edger_method, 0 for quasi-likelihood (QL) F-test, 1 for likelihood ratio test, 2 for exact test
#@param list edger_contrast, results for which groups comparison

edger_identify_degs <- function() {
	#convert to integer
	read_counts <- round(read_counts)

	#get comparison group information
	info <- sample_info[colnames(read_counts),]
	groups <- info[, edger_group]

	#make deg list with group
	y <- DGEList(counts=read_counts, group=groups)

	#filter out lowly expressed genes
	keep <- filterByExpr(y)
	y <- y[keep, , keep.lib.sizes=FALSE]

	#perform the TMM normalization
	y <- normLibSizes(y)

	#make edgeR design matrix
	edger_design <<- model.matrix(as.formula(edger_design), data=info)
	colnames(edger_design) <- levels(y$samples$group)

	#estimate common dispersion and tagwise dispersions
	edger_count <<- estimateDisp(y, edger_design)
}

edger_make_contrast <- function() {
	cmp <- paste(edger_contrast, collapse='-')
	con <- makeContrasts(cmp, levels=edger_design)
	return(con)
}

edger_plot_degs <- function() {
	tags <- sig_degs[,1]
	plotSmear(edger_degs, de.tags=tags)
}

deger_sig_degs <- function() {
	sig <- topTags(et, n=Inf)
	sig_degs <<- sig$table[sig$table$FDR < edger_fdr & sig$table$logFC>=edger_logfc,]
}

edger_quasi_likelihood <- function() {
	fit <- glmQLFit(edger_count, edger_design)
	con <- edger_make_contrast()
	edger_degs <<- glmQLTest(fit, contrast=con)
}

edger_likelihood_ratio <- function() {
	fit <- glmFit(edger_count, edger_design)
	con <- edger_make_contrast()
	edger_degs <<- glmLRT(fit, contrast=con)
}

edger_exact_test <- function() {
	edger_degs <<- exactTest(edger_count, pair=edger_contrast)
}

edger_analysis_pipeline <- function() {

}

edger_show_degs <- function() {
	
}
