# Refresh data to use in test_admrs.R

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

admrs_output_mean_max <- admrs(wt,  mut_dl_sigma = 0.3, mut_sd_min = 1,  option = "mean_max")
usethis::use_data(admrs_output_mean_max,  overwrite = TRUE)

admrs_output_max_max <- admrs(wt,  mut_dl_sigma = 0.3, mut_sd_min = 1,  option = "max_max")
usethis::use_data(admrs_output_max_max,  overwrite = TRUE)



