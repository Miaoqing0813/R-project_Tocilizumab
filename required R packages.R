
install.packages(c("devtools", "dplyr", "ggplot2", "reshape2", "data.table", 
                   "hesim", "here", "flexsurv", "IPDfromKM", "mixcure", "zoo"
                   "kableExtra", "doRNG", "doParallel", "survminer", "tm"))

library(devtools)
devtools::install_github("Pechli-Lab/SurvdigitizeR")
devtools::install_github("DARTH-git/darthtools", dependencies = TRUE)
remotes::install_github("nx10/httpgd")

library(dplyr)
library(knitr)
library(ggplot2)
library(data.table)
library (here)
library(flexsurv)
library(SurvdigitizeR)
library(IPDfromKM)
library(survival)
library(mixcure)
library(httpgd)
library(kableExtra)
library (future)
library(doRNG)
library(doParallel)
library(darthtools)
library(survminer)
library(tm)
library(ggsurvfit)
library(zoo)

