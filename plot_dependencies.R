library(here)
library(mvbutils)

remove(list = ls())
# source("R/enm.R")
# source("R/penm.R")
# source("R/penm_analysis.R")

source("R/prs_new.R")
source("R/prs_fast.R")

fw <- foodweb()

foodweb(prune = c("prs_all.new"))
foodweb(prune = c("prs_all.fast"))



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


