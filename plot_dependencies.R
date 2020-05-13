library(here)
library(mvbutils)

remove(list = ls())
source("R/enm.R")
source("R/penm.R")

fw <- foodweb()

foodweb(prune =c("set_enm"))

foodweb(prune =c("get_mutant_site"))

foodweb(prune =c("get_mutant_site_lfenm"))

foodweb(prune =c("get_mutant_site_sclfenm"),
        generics = c("dij_edge", "sdij_edge"))

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


