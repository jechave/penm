library(here)
library(tidyverse)
library(bio3d)


load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "wt.rda"))
load(test_path("fixtures", "mut_lf.rda"))
load(test_path("fixtures", "mut_qf.rda"))


test_that("set_enm gets wt ok", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE), wt)
})

test_that("get_mutant_site gets mut_lf", {
  expect_equal(
    get_mutant_site(wt, site_mut = 80, mutation = 1,
                    mut_model = "lfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_lf)
})

test_that("get_mutant_site gets mut_qf", {
  skip("Skip sclfenm test until sclefnm is fixed")
  expect_equal(
    get_mutant_site(wt,  site_mut = 80, mutation = 1,
                    mut_model = "sclfenm", mut_sd_min = 1, mut_dl_sigma = 0.3),
    mut_qf)
})
