# From the R console, run `source("install.packages.for.BI347.R")`

{
  start.time <- Sys.time()

  # Begin the package installations
  install.packages('stringr')
  install.packages('magrittr')
  install.packages('Momocs')
  install.packages('devtools')
  devtools::install_github("mlcollyer/RRPP")

  cat("\nDone!\n\n")

  # Check that everything was installed
  x <- c("stringr","magrittr","Momocs","devtools","RRPP")
  if (all(x %in% as.data.frame(installed.packages())[,1])) {
    cat("All packages successfully install.\n")
  } else {
    i <- which(!(x %in% as.data.frame(installed.packages())[,1]) )
    cat("ATTENTION - The following packages were not installed: \n")
    cat(paste0(x[i],collapse = "\n"),"\n\n")
  }
  print(Sys.time() - start.time)
}
