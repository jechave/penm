library(here)
library(mvbutils)

remove(list = ls())
source("R/enm.R")

source("R/penm.R")
source("R/penm_analysis.R")

source("R/prs.R")
source("R/prs_fast.R")

fw <- foodweb()

foodweb(prune = c("calculate_dfij.prs"))
foodweb(prune = c("calculate_dfij.fast"))
foodweb(prune = c("calculate_dfj.prs"))
foodweb(prune = c("calculate_dfj.fast"))


foodweb(prune = c("calculate_dfj.prs", "calculate_dfj.fast"))
foodweb(prune = c("calculate_dfj", "calculate_dfj.fast"))
foodweb(prune = c("calculate_dfij", "calculate_dfij.fast"))
foodweb(prune = c("calculate_dfnj", "calculate_dfnj.fast"))

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


