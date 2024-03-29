## ----------------------------------------------------------------------------------------------------------------------
# Test analytic Double Mutational Response Scanning function `admrs`

# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)
library(tictoc)
library(Matrix)

load(here("data/wt.rda"))
load(here("data/admrs_output_mean_max.rda"))
load(here("data/admrs_output_max_max.rda"))

# Test analytic dmrs

admrs_output_mean_max_test <- admrs(wt,  mut_dl_sigma = 0.3, mut_sd_min = 1,  option = "mean_max")

test_that("admrs with option mean_max is ok", {
  expect_equal(admrs_output_mean_max_test, admrs_output_mean_max)
  })


admrs_output_max_max_test <- admrs(wt,  mut_dl_sigma = 0.3, mut_sd_min = 1,  option = "max_max")

test_that("admrs with option max_max  is ok", {
  expect_equal(admrs_output_max_max_test, admrs_output_max_max)
  })

# Test sdmrs (simulation-based dmrs)

load(here("data/sdmrs_output_mean_max.rda"))
load(here("data/sdmrs_output_max_max.rda"))


sdmrs_output_mean_max_test <- sdmrs(wt, nmut = 5, mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234, option = "mean_max")

test_that("sdmrs with option mean_max is ok", {
  expect_equal(sdmrs_output_mean_max_test$dmrs_matrix, sdmrs_output_mean_max$dmrs_matrix)
})



sdmrs_output_max_max_test <- sdmrs(wt, nmut = 5, mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1234, option = "max_max")

test_that("sdmrs with option max_max  is ok", {
  expect_equal(sdmrs_output_max_max_test$dmrs_matrix, sdmrs_output_max_max$dmrs_matrix)
})




