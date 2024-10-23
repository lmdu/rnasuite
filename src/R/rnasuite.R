RNASUITE.PLOTS <- list()

RnasuiteGetPlot <- function(id) {
	id <- as.character(id)
	plot <- RNASUITE.PLOTS[[ id ]]
	return(plot)
}

RnasuiteGetPlotIdByName <- function(name) {
	res <- Filter(function(x) x$name == name, RNASUITE.PLOTS)
	
	if (length(res) > 0) {
		id <- res[[1]]$id
	} else {
		id <- NULL
	}

	return(id)
}

RnasuiteGetPlotCurrentId <- function() {
	return(unigd::ugd_id()$id)
}

RnasuiteDeletePlot <- function(id) {
	id <- as.character(id)
	RNASUITE.PLOTS[[ id ]] <<- NULL
}

RnasuiteSavePlot <- function(old, new, name, plot, data) {
	if (!is.null(old)) {
		RnasuiteDeletePlot(old)
	}

	index <- as.character(new)
	RNASUITE.PLOTS[[ index ]] <<- list(
		id = new,
		name = name,
		plot = plot,
		data = data
	)
}

#convert between pandas and dataframe
RnasuitePandasToDataframe <- function(data) {
	return(py_to_r(r_to_py(data)))
}

RnasuiteDataframeToPandas <- function(data) {
	return(r_to_py(data))
}
