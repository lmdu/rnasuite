#@param dataframe read_counts, raw read count expression data
#@param dataframe sample_info, sample meta and group information
#@param str edger_design, design formula for deg identification
#@param float edger_fdr, FDR value for deg identification
#@param float edger_logfc, log2foldchange for deg identification
#@param str edger_group, find DEGs between given group
#@param int edger_replicate, 0 for with biological replicates, 1 for with on replicates
#@param int edger_method, 0 for quasi-likelihood (QL) F-test, 1 for likelihood ratio test, 2 for exact test
#@param list edger_contrast, results for which groups comparison
#@param float edger_bcv, bcv value for dispersion
library(edgeR)

edger_data_prepare <- function() {
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
	edger_model <<- model.matrix(as.formula(edger_design), data=info)
	#colnames(edger_design) <- levels(y$samples$group)
	colnames(edger_model) <- gsub(paste0('^', edger_group), '', colnames(edger_model))

	#estimate common dispersion and tagwise dispersions
	if (edger_replicate == 0) {
		if (edger_method < 2) {
			edger_count <<- estimateDisp(y, edger_model)
		} else {
			edger_count <<- estimateDisp(y)
		}
	} else {
		edger_count <<- y
	}
}

edger_make_contrast <- function() {
	cmp <- paste(edger_contrast, collapse='-')
	con <- makeContrasts(cmp, levels=edger_model)
	return(con)
}

edger_with_replicates <- function() {
	if (edger_method < 2) {
		con <- edger_make_contrast()

		if (edger_method == 0) {
			fit <- glmQLFit(edger_count, edger_model)
			edger_degs <<- glmQLTest(fit, contrast=con)
		} else {
			fit <- glmFit(edger_count, edger_model)
			edger_degs <<- glmLRT(fit, contrast=con)
		}
	} else {
		edger_degs <<- exactTest(edger_count, pair=edger_contrast)
	}
}

edger_without_replicates <- function() {
	if (edger_method == 1) {
		fit <- glmFit(edger_count, edger_model, dispersion=edger_bcv^2)
		con <- edger_make_contrast()
		edger_degs <<- glmLRT(fit, contrast=con)
	} else {
		edger_degs <<- exactTest(edger_count, pair=edger_contrast, dispersion=edger_bcv^2)
	}
}

edger_identify_degs <- function() {
	if (edger_replicate == 0) {
		edger_with_replicates()
	} else {
		edger_without_replicates()
	}
}

edger_sig_degs <- function() {
	sig <- topTags(et, n=Inf)
	sig_degs <<- sig$table[sig$table$FDR < edger_fdr & sig$table$logFC>=edger_logfc,]
}

edger_plot_degs <- function() {
	tags <- sig_degs[,1]
	plotSmear(edger_degs, de.tags=tags)
}

edger_return_degs <- function() {
	out <- list(
		normal_count = cpm(edger_count),
		degs_list = sig_degs,
		degs_versus = deseq_contrast
	)
	return(out)
}

edger_analysis_pipeline <- function() {
	edger_data_prepare()
	edger_identify_degs()
	edger_sig_degs()
	edger_plot_degs()
	edger_return_degs()
}

edger_show_degs <- function() {

}
