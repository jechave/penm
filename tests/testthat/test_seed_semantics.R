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

test_that("same (site, mutation) and same ensemble give the identical mutation", {
  a <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)
  b <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)

  expect_identical(dlij(a), dlij(b))
  # Not vacuous: something was actually perturbed, so the line above is not
  # comparing two all-zero vectors.
  expect_true(any(dlij(a) != 0))
})

test_that("the identical mutation survives unrelated RNG use in between", {
  # A mutant's identity must not depend on session history.
  a <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)
  set.seed(99)
  rnorm(1000)
  b <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)

  expect_identical(dlij(a), dlij(b))
})

test_that("a different ensemble is a different realization", {
  a <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)
  b <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 2L)

  # The perturbed contacts are a property of the site, not of the realization,
  # so the edge set must NOT move -- only the magnitudes.
  expect_identical(which(dlij(a) != 0), which(dlij(b) != 0))
  expect_false(isTRUE(all.equal(dlij(a), dlij(b))))
})

test_that("different (site, mutation) in one ensemble are different mutations", {
  m3 <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)
  m4 <- get_mutant_site(wt, site_mut = 80, mutation = 4, mut_sd_min = 1, ensemble = 1L)
  expect_false(isTRUE(all.equal(dlij(m3), dlij(m4))))

  # Comparing WHICH edges each site perturbs would prove nothing: that is set
  # by the contact graph, so sites 80 and 81 differ however broken the key is.
  # The key controls the magnitudes, so compare those -- on one site, holding
  # the graph fixed and varying only site_mut's contribution to the key.
  s80 <- get_mutant_site(wt, site_mut = 80, mutation = 3, mut_sd_min = 1, ensemble = 1L)
  expect_false(isTRUE(all.equal(
    penm:::mut_seed(1L, 80, 3),
    penm:::mut_seed(1L, 81, 3)
  )))
  # and the draws that follow from those keys differ
  draw <- function(k) { set.seed(k); rnorm(5) }
  expect_false(isTRUE(all.equal(draw(penm:::mut_seed(1L, 80, 3)),
                                draw(penm:::mut_seed(1L, 81, 3)))))
})

test_that("a malformed ensemble is rejected, not turned into a realization", {
  # Each of these used to produce a valid-looking mutant: the key is built with
  # paste(), which stringifies anything. NULL gave the key "-80-3".
  expect_error(get_mutant_site(wt, 80, 3, ensemble = NULL),   "single non-missing integer")
  expect_error(get_mutant_site(wt, 80, 3, ensemble = NA),     "single non-missing integer")
  expect_error(get_mutant_site(wt, 80, 3, ensemble = "banana"), "single non-missing integer")
  expect_error(get_mutant_site(wt, 80, 3, ensemble = c(1, 2)), "single non-missing integer")
  expect_error(get_mutant_site(wt, 80, 3, ensemble = 1.5),    "single non-missing integer")

  # Fires on the mutation = 0 path too, which returns early and never reaches
  # mut_seed(). Accepting a bad value there while erroring at mutation = 1
  # would be a trap.
  expect_error(get_mutant_site(wt, 80, 0, ensemble = NULL),   "single non-missing integer")

  # Valid values, integer and double, are accepted.
  expect_s3_class(get_mutant_site(wt, 80, 3, mut_sd_min = 1, ensemble = 7L), "prot")
  expect_s3_class(get_mutant_site(wt, 80, 3, mut_sd_min = 1, ensemble = 7),  "prot")
})

test_that("mut_seed rejects a malformed ensemble when called directly", {
  # Not redundant with the boundary check: penmscan and any direct penm:::
  # caller reach mut_seed() without passing through get_mutant_site().
  expect_error(penm:::mut_seed(NULL, 80, 3), "single non-missing integer")
  expect_error(penm:::mut_seed(NA,   80, 3), "single non-missing integer")
})

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
