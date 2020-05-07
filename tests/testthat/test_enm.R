library(here)

load(here("data/pdb_2acy_A.rda"))
load(here("data/prot_2acy_A_ming_wall_ca.rda"))


test_that("set_enm gets prot  equal to prot_2acy_A", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE),
               prot_2acy_A_ming_wall_ca)
})
