library(here)
library(mvbutils)

remove(list = ls())
# source("R/enm.R")
# source("R/penm.R")
# source("R/penm_analysis.R")

list.files(path = "R")
for (f in list.files(path = "R")) {
  source(here("R",f))
}

result <- foodweb(plotting = FALSE)

# The following line returns the number of times each function is called
res <- sapply(rownames(result$funmat), function(n) length(callers.of(n, result)))


# Get those functions that are never called:
names(res[res==0])


# callers.of("enm_v_min", fw)
# callers.of("calculate_dvm", fw)
#
# callers.of("enm_v_min", fw)
# callers.of("get_mutant_site", fw)
#
# foodweb(prune = c("prs_all.new"))
# foodweb(prune = c("prs_all.fast"))




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


