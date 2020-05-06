library(mvbutils)

remove(list = ls())
source("R/enm.R")
source("R/enm_prot_getters.R")

fw <- foodweb()

callers.of("set_enm", fw)
callees.of("set_enm", fw)
