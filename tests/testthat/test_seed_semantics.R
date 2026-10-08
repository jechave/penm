# Semantics of the mutant-identity argument: ensemble, set by set_enm() and
# used by get_mutant_site().
#
# test_seed.R pins the properties of the hash key itself. This file pins that
# those properties reach the perturbations a user actually gets: that
# (ensemble, site_mut, mutation) names one reproducible mutation.
#
# These go through get_mutant_site(), not mut_seed(), and compare graph$lij --
# the quantity a mutation perturbs. Comparing whole prot objects would pass or
# fail for reasons unrelated to the seeding.

library(penm)

load(test_path("fixtures", "pdb_2acy_A.rda"))

wt_in_ensemble <- function(ensemble) {
  set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
          mut_sd_min = 1, ensemble = ensemble)
}
wt <- wt_in_ensemble(1L)
wt_2 <- wt_in_ensemble(2L)

dlij <- function(mut) mut$graph$lij - wt$graph$lij

test_that("same (site, mutation) and same ensemble give the identical mutation", {
  a <- get_mutant_site(wt, site_mut = 80, mutation = 3)
  b <- get_mutant_site(wt, site_mut = 80, mutation = 3)

  expect_identical(dlij(a), dlij(b))
  # Not vacuous: something was actually perturbed, so the line above is not
  # comparing two all-zero vectors.
  expect_true(any(dlij(a) != 0))
})

test_that("the identical mutation survives unrelated RNG use in between", {
  # A mutant's identity must not depend on session history.
  a <- get_mutant_site(wt, site_mut = 80, mutation = 3)
  set.seed(99)
  rnorm(1000)
  b <- get_mutant_site(wt, site_mut = 80, mutation = 3)

  expect_identical(dlij(a), dlij(b))
})

test_that("a different ensemble is a different realization", {
  a <- get_mutant_site(wt, site_mut = 80, mutation = 3)
  b <- get_mutant_site(wt_2, site_mut = 80, mutation = 3)

  # The perturbed contacts are a property of the site, not of the realization,
  # so the edge set must NOT move -- only the magnitudes.
  expect_identical(which(dlij(a) != 0), which(dlij(b) != 0))
  expect_false(isTRUE(all.equal(dlij(a), dlij(b))))
})

test_that("different (site, mutation) in one ensemble are different mutations", {
  m3 <- get_mutant_site(wt, site_mut = 80, mutation = 3)
  m4 <- get_mutant_site(wt, site_mut = 80, mutation = 4)
  expect_false(isTRUE(all.equal(dlij(m3), dlij(m4))))

  # Comparing WHICH edges each site perturbs would prove nothing: that is set
  # by the contact graph, so sites 80 and 81 differ however broken the key is.
  # The key controls the magnitudes, so compare those -- on one site, holding
  # the graph fixed and varying only site_mut's contribution to the key.
  s80 <- get_mutant_site(wt, site_mut = 80, mutation = 3)
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
  expect_error(wt_in_ensemble(NULL),     "single non-missing integer")
  expect_error(wt_in_ensemble(NA),       "single non-missing integer")
  expect_error(wt_in_ensemble("banana"), "single non-missing integer")
  expect_error(wt_in_ensemble(c(1, 2)),  "single non-missing integer")
  expect_error(wt_in_ensemble(1.5),      "single non-missing integer")

  # Valid values, integer and double, are accepted, and name the same realization.
  m_integer <- get_mutant_site(wt_in_ensemble(7L), 80, 3)
  m_double  <- get_mutant_site(wt_in_ensemble(7),  80, 3)
  expect_s3_class(m_double, "prot")
  expect_identical(m_double$graph$lij, m_integer$graph$lij)
  expect_identical(get_enm_param(m_double)$ensemble, 7L)
})

test_that("a mutant inherits the ensemble, so a trajectory stays in one realization", {
  # mutating a mutant: the second mutation is the one wt would get at that site
  m1 <- get_mutant_site(wt, site_mut = 80, mutation = 3)
  m2 <- get_mutant_site(m1, site_mut = 20, mutation = 5)
  direct <- get_mutant_site(wt, site_mut = 20, mutation = 5)
  expect_identical(get_enm_param(m2)$ensemble, 1L)
  # equal, not identical: (lij + d) - lij need not be d to the last bit
  expect_equal(m2$graph$lij - m1$graph$lij, dlij(direct))
  expect_true(any(dlij(direct) != 0))
})

test_that("a wt built by an earlier penm, without mut_model, is rejected", {
  old_wt <- wt
  old_wt$param <- old_wt$param[c("node", "model", "d_max")]
  expect_error(get_mutant_site(old_wt, 80, 3), "Rebuild it with set_enm")
})

test_that("mut_seed rejects a malformed ensemble when called directly", {
  # Not redundant with the boundary check: a direct penm::: caller reaches
  # mut_seed() without passing through set_enm().
  expect_error(penm:::mut_seed(NULL, 80, 3), "single non-missing integer")
  expect_error(penm:::mut_seed(NA,   80, 3), "single non-missing integer")
})

test_that("get_mutant_site does not disturb the caller's RNG stream", {
  # set.seed() writes .Random.seed in the global environment, so seeding a
  # mutant's perturbations also silently reseeds whatever the caller was doing.
  set.seed(42)
  expected <- rnorm(3)

  set.seed(42)
  invisible(get_mutant_site(wt, site_mut = 80, mutation = 1))

  expect_equal(rnorm(3), expected)
})
