send_val_to_r <- function(var, data) {
	assign(var, data, envir=.GlobalEnv)
}

send_df_to_r <- function(var, data) {
	assign(var, py_to_r(r_to_py(data)), envir=.GlobalEnv)
}
