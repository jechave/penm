# Refresh data to use in test_penm.R

## ----------------------------------------------------------------------------------------------------------------------
# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)


## ----------------------------------------------------------------------------------------------------------------------
load(here("data/pdb_2acy_A.rda"))

wt_sc <- set_enm(pdb_2acy_A, node = "sc", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
usethis::use_data(wt_sc,  overwrite = TRUE)

## ----------------------------------------------------------------------------------------------------------------------

mut_sc_lf  <- get_mutant_site(wt, site_mut = 80, mutation = 1,
                           mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1)

usethis::use_data(mut_sc_lf,  overwrite = TRUE)

## ----------------------------------------------------------------------------------------------------------------------
skip <-  TRUE
if (!skip) {
  mut_sc_qf <- get_mutant_site(wt, site_mut = 80, mutation = 1,
                            mut_model = "sclfenm", mut_dl_sigma = 0.3, mut_sd_min = 1)
  usethis::use_data(mut_sc_qf, overwrite = TRUE)
}






