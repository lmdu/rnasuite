RNASUITE_PLOTS <- list()
#RNASUITE_DATAS <- list()
#RNASUITE_NAMES <- list()

rnasuite_get_plot <- function(id) {
	if (is.numeric(id)) {
		id <- as.character(id)
		plot <- RNASUITE_PLOTS[[ id ]]
	} else {
		res <- Filter(function(x) x$name == id, RNASUITE_PLOTS)
		if (length(res) > 0) {
			plot <- (res[[1]])
		} else {
			plot <- NULL
		}
	}

	return(plot)
}

rnasuite_del_plot <- function(id) {
	id <- as.character(id)
	RNASUITE_PLOTS[[ id ]] <<- NULL
}

rnasuite_put_plot <- function(old, new, name, plot, data) {
	if (!is.null(old)) {
		rnasuite_del_plot(old)
	}

	index <- as.character(new)
	RNASUITE_PLOTS[[ index ]] <<- list(
		id = new,
		name = name,
		plot = plot,
		data = data
	)
}

rnasuite_get_plot_id <- function() {
	return(unigd::ugd_id()$id)
}

rnasuite_pandas_to_dataframe <- function(data) {
	return(py_to_r(r_to_py(data)))
}

rnasuite_dataframe_to_pandas <- function(data) {
	return(r_to_py(data))
}
