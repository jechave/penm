library(mvbutils)

remove(list = ls())
source("plot_enm.R")
source("R/enm_analysis.R")

fw <- foodweb()

callers.of("matrix_to_tibble", fw)
callees.of("matrix_to_tibble", fw)
