# From the R console, run `source("install.packages.for.BI377.R`

{
  start.time <- Sys.time()
  
  # Begin the package installations
  install.packages('tidyverse')
  install.packages('magrittr')
  install.packages('ggpubr')
  install.packages('ggrepel')
  install.packages('ggplotify')
  install.packages('jpeg') # Must install before geomorph
  install.packages('viridis')
  install.packages('qrcode')
  install.packages('beeswarm')
  install.packages('car')
  install.packages('vegan')
  install.packages('dunn.test')
  install.packages('factoextra')
  install.packages('nlme')
  install.packages('lme4')
  # install.packages('ape')
  install.packages('phytools')
  # install.packages('phylotools') # Recently kicked off CRAN!
  install.packages('palmerpenguins')
  install.packages('geomorph')
  install.packages('Momocs')
  
  # install.packages('curl')
  # install.packages('usethis')
  install.packages('devtools')
  devtools::install_github("mlcollyer/RRPP")
  devtools::install_github("helixcn/phylotools", build_vignettes = TRUE)
  devtools::install_github("aphanotus/borealis", build_vignettes = TRUE)

  cat("\nDone!\n\n")
  
  # Check that everything was installed
  x <- c("dplyr","ggplot2","stringr","magrittr","ggpubr","ggrepel","ggplotify","viridis",
         "qrcode","beeswarm","car","vegan","dunn.test","factoextra","nlme","lme4",
         "ape","phytools","phylotools","palmerpenguins","geomorph","RRPP","Momocs",
         "devtools","borealis")
  if (all(x %in% as.data.frame(installed.packages())[,1])) {
    cat("All packages successfully install.\n")
  } else {
    i <- which(!(x %in% as.data.frame(installed.packages())[,1]) )
    cat("ATTENTION - The following packages were not installed: \n")
    cat(paste0(x[i],collapse = "\n"),"\n\n")
  }
  print(Sys.time() - start.time)
}
