library(here)
library(mvbutils)

remove(list = ls())
source("R/enm.R")
source("R/enm_energy.R")
source("R/enm_analysis.R")
source("R/enm_kij_functions.R")
source("R/enm_prot_getters.R")
source("R/penm.R")
source("R/penm_analysis.R")
source("R/sdmrs.R")
source("R/admrs.R")

source("R/mrs.R")
source("R/smrs.R")
source("R/amrs.R")

fw <- foodweb()

foodweb(prune = c("mrs_all"))

foodweb(prune = c("smrs_all"))

foodweb(prune = c("amrs_all"))

foodweb(prune = c("sdmrs"))


# callers.of("get_mutant_site", fw)
# callees.of("get_mutant_site", fw)
# callees.of(funs = callees.of("get_mutant_site", fw), fw)
#
# library(tidyverse)
#
# "get_mutant_site" %>%
#   callees.of(fw = fw) %>%
#   callees.of(fw = fw) %>%
#   callees.of(fw = fw)


