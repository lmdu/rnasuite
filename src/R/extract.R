rnasuite_generic_extract_degs <- function(treatment, control, ...) {
	tool <- RNASUITE_DEGTOOL

	if (tool == 'deseq') {
		out <- rnasuite_deseq_extract_degs(treatment, control)
	}

	return(out)
}
