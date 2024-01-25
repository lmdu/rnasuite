RNASUITE_PLOTS <- list()
RNASUITE_DATAS <- list()
RNASUITE_NAMES <- list()

rnasuite_get_plot <- function(id) {
	id <- as.character(id)
	return(RNASUITE_PLOTS[[ id ]])
}

rnasuite_get_data <- function(id) {
	id <- as.character(id)
	return(RNASUITE_DATAS[[ id ]])
}

rnasuite_get_name <- function(id) {
	id <- as.character(id)
	return(RNASUITE_NAMES[[ id ]])
}

rnasuite_get_id <- function(name) {
	id = names(which(RNASUITE_NAMES == name))

	if (length(id) > 0) {
		id = as.integer(id)
	} else {
		id = NULL
	}

	return(id)
}

rnasuite_put_plot <- function(old, new, name, plot, data) {
	if (!is.null(old)) {
		old <- as.character(old)
		RNASUITE_PLOTS[[ old ]] <<- NULL
		RNASUITE_DATAS[[ old ]] <<- NULL
		RNASUITE_NAMES[[ old ]] <<- NULL
	}

	new <- as.character(new)
	RNASUITE_NAMES[[ new ]] <<- name
	RNASUITE_PLOTS[[ new ]] <<- plot
	RNASUITE_DATAS[[ new ]] <<- data
}

#make plot id start from one
plot.new()
