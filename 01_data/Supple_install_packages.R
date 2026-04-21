
install.packages(c("devtools", "dplyr", "ggplot2", "reshape2", "data.table", 
                   "hesim", "here", "flexsurv", "IPDfromKM", "mixcure", "zoo",
                   "kableExtra", "doRNG", "doParallel", "survminer", "tm", "survfit"))
install.packages("openxlsx", dependencies = TRUE)

library(devtools)
devtools::install_github("Pechli-Lab/SurvdigitizeR")
devtools::install_github("DARTH-git/darthtools", dependencies = TRUE)
remotes::install_github("nx10/httpgd")




