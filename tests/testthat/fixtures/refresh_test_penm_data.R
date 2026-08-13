# Refresh data to use in test_penm.R

## ----------------------------------------------------------------------------------------------------------------------
# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)


## ----------------------------------------------------------------------------------------------------------------------
load(here("tests/testthat/fixtures/pdb_2acy_A.rda"))

wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
save(wt, file = here("tests/testthat/fixtures/wt.rda"))

## ----------------------------------------------------------------------------------------------------------------------

mut_lf  <- get_mutant_site(wt, site_mut = 80, mutation = 1,
                           mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1)

save(mut_lf, file = here("tests/testthat/fixtures/mut_lf.rda"))

## ----------------------------------------------------------------------------------------------------------------------
skip <-  TRUE
if (!skip) {
  mut_qf <- get_mutant_site(wt, site_mut = 80, mutation = 1,
                            mut_model = "sclfenm", mut_dl_sigma = 0.3, mut_sd_min = 1)
  save(mut_qf, file = here("tests/testthat/fixtures/mut_qf.rda"))
}






