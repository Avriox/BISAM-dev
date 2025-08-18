library(Rcpp)
setwd("~/Documents/Uni/WU/BISAM")
Rcpp::compileAttributes(verbose = TRUE)
install.packages("./", repos = NULL, type = "source", verbose = TRUE, clean=TRUE)
