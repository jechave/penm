load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "wt.rda"))
load(test_path("fixtures", "mut_lf.rda"))


test_that("set_enm gets wt ok", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5), wt)
})

test_that("get_mutant_site gets mut_lf", {
  expect_equal(
    get_mutant_site(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                            mut_dl_sigma = 0.3, mut_sd_min = 1),
                    site_mut = 80, mutation = 1),
    mut_lf)
})
