.First.lib <- function(lib,pkg) {
	library.dynam("pscn",pkg,lib)
   	cat("package pscn loaded\n")
}
