library(here)

load(here("data/pdb_2acy_A.rda"))
load(here("data/prot_2acy_A.rda"))


test_that("set_prot gets prot  equal to prot_2acy_A", {
  expect_equal(set_prot(pdb_2acy_A, node = "sc", model = "ming_wall", d_max = 10.5, frustrated = FALSE),
               prot_2acy_A)
})
