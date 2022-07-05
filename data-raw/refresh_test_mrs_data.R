# Refresh data to use in test_amrs.R

## ----------------------------------------------------------------------------------------------------------------------
# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)
library(tictoc)
library(Matrix)


## ----------------------------------------------------------------------------------------------------------------------
load(here("data/wt.rda"))

mrs_all_output <- mrs_all(wt, nmut_per_site = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
usethis::use_data(mrs_all_output,  overwrite = TRUE)

smrs_all_output <- smrs_all(wt, nmut_per_site = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
usethis::use_data(smrs_all_output,  overwrite = TRUE)

amrs_all_output <- amrs_all(wt, mut_dl_sigma = 0.3, mut_sd_min = 1)
usethis::use_data(amrs_all_output,  overwrite = TRUE)




