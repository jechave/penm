# Semantics of the mutant-identity argument at the get_mutant_site() boundary.
#
# test_seed.R pins the properties of the hash key itself. This file pins that
# those properties reach the perturbations a user actually gets: that
# (ensemble, site_mut, mutation) names one reproducible mutation.
#
# These go through get_mutant_site(), not mut_seed(), and compare graph$lij --
# the quantity a mutation perturbs. Comparing whole prot objects would pass or
# fail for reasons unrelated to the seeding.

library(penm)

load(test_path("fixtures", "wt.rda"))

dlij <- function(mut) mut$graph$lij - wt$graph$lij

test_that("get_mutant_site does not disturb the caller's RNG stream", {
  # set.seed() writes .Random.seed in the global environment, so seeding a
  # mutant's perturbations also silently reseeds whatever the caller was doing.
  set.seed(42)
  expected <- rnorm(3)

  set.seed(42)
  invisible(get_mutant_site(wt, site_mut = 80, mutation = 1,
                            mut_model = "lfenm", mut_sd_min = 1))

  expect_equal(rnorm(3), expected)
})
