RNASUITE_PLOTS <- list()
RNASUITE_DATAS <- list()
RNASUITE_NAMES <- list()

rnasuite_get_plot <- function(id) {
	return(RNASUITE_PLOTS[[ id ]])
}

rnasuite_get_data <- function(id) {
	return(RNASUITE_DATAS[[ id ]])
}

rnasuite_get_name <- function(id) {
	return(RNASUITE_NAMES[[ id ]])
}

rnasuite_put_plot <- function(old, new, name, plot, data=NULL) {
	if (!is.null(old)) {
		RNASUITE_PLOTS[[ old ]] <<- NULL
		RNASUITE_DATAS[[ old ]] <<- NULL
		RNASUITE_NAMES[[ old ]] <<- NULL
	}

	RNASUITE_NAMES[[ new ]] <<- name
	RNASUITE_PLOTS[[ new ]] <<- plot

	if (!is.null(data)) {
		RNASUITE_DATAS[[ new ]] <<- data
	}
}



