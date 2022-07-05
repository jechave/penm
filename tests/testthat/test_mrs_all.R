## ----------------------------------------------------------------------------------------------------------------------
# Test analytic Double Mutational Response Scanning function `amrs`

# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)
library(tictoc)
library(Matrix)

load(here("data/wt.rda"))

load(here("data/mrs_all_output.rda"))
mrs_all_output_test <- mrs_all(wt, nmut_per_site = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
test_that("mrs with option mean_max is ok", {
  expect_equal(mrs_all_output_test, mrs_all_output)
  })

load(here("data/smrs_all_output.rda"))
smrs_all_output_test <- smrs_all(wt, nmut_per_site = 5, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234)
test_that("smrs with option mean_max is ok", {
  expect_equal(smrs_all_output_test, smrs_all_output)
})

load(here("data/amrs_all_output.rda"))
amrs_all_output_test <- amrs_all(wt, mut_dl_sigma = 0.3, mut_sd_min = 1)
test_that("amrs with option mean_max is ok", {
  expect_equal(amrs_all_output_test, amrs_all_output)
  })



