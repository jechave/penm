library(here)

load(here("data/pdb_2acy_A.rda"))
load(here("data/wt.rda"))
load(here("data/mut_lf.rda"))
load(here("data/mut_qf.rda"))


test_that("set_enm gets wt ok", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE), wt)
})

test_that("get_mutant_site gets mut_lf", {
  expect_equal(
    get_mutant_site(wt, site_mut = 80, mutation = 1,
                    wt0 = wt,
                    mut_sd_min = 1, dl_sigma = 0.3, update_enm = FALSE,
                    model = "ming_wall", d_max = 10.5, frustrated = FALSE),
    mut_lf)
})

test_that("get_mutant_site gets mut_qf", {
  expect_equal(
    get_mutant_site(wt, site_mut = 80, mutation = 1,
                    wt0 = wt,
                    mut_sd_min = 1, dl_sigma = 0.3, update_enm = TRUE,
                    model = "ming_wall", d_max = 10.5, frustrated = FALSE),
    mut_qf)
})
