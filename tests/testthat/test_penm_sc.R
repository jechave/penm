library(here)
library(tidyverse)
library(jefuns)
library(bio3d)


load(here("data/pdb_2acy_A.rda"))
load(here("data/wt_sc.rda"))
load(here("data/mut_sc_lf.rda"))
# load(here("data/mut_sc_qf.rda"))


test_that("set_enm gets wt_sc ok", {
  expect_equal(set_enm(pdb_2acy_A, node = "sc", model = "ming_wall", d_max = 10.5, frustrated = FALSE),
               wt_sc)
})

test_that("get_mutant_site gets mut_sc_lf", {
  expect_equal(
    get_mutant_site(wt_sc, site_mut = 80, mutation = 1,
                    mut_model = "lfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_sc_lf)
})

test_that("get_mutant_site gets mut_qf", {
  skip("Skip sclfenm test until sclefnm is fixed")
  expect_equal(
    get_mutant_site(wt_sc,  site_mut = 80, mutation = 1,
                    mut_model = "sclfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_sc_qf)
})
