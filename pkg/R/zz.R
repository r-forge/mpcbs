.First.lib <- function(lib,pkg) {
   library.dynam("msscan",pkg,lib)
   library(fields)
   cat("roots 0.1-1 loaded\n")
}

