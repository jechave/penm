library(penm)

load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "prot_2acy_A_ming_wall_ca.rda"))


test_that("set_enm gets prot  equal to prot_2acy_A", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5),
               prot_2acy_A_ming_wall_ca)
})

test_that("set_enm works with beta carbon nodes", {
  # Test that CB nodes can be created
  prot_cb <- set_enm(pdb_2acy_A, node = "cb", model = "anm", d_max = 12.0)

  # Check that the object is created properly
  expect_s3_class(prot_cb, "prot")
  expect_equal(prot_cb$param$node, "cb")

  # Check that we have the right number of nodes
  expect_equal(prot_cb$nodes$nsites, length(prot_cb$nodes$pdb_site))

  # Check that coordinates are defined (not all NA)
  expect_false(all(is.na(prot_cb$nodes$xyz)))

  # Check that the ENM components are created
  expect_false(is.null(prot_cb$kmat))
  expect_false(is.null(prot_cb$nma))
  expect_false(is.null(prot_cb$graph))
})

test_that("set_enm accepts both 'cb' and 'beta' for beta carbon nodes", {
  # Test that both 'cb' and 'beta' work
  prot_cb1 <- set_enm(pdb_2acy_A, node = "cb", model = "anm", d_max = 12.0)
  prot_cb2 <- set_enm(pdb_2acy_A, node = "beta", model = "anm", d_max = 12.0)

  # Both should produce the same result
  expect_equal(prot_cb1$nodes$xyz, prot_cb2$nodes$xyz)
  expect_equal(prot_cb1$nodes$bfactor, prot_cb2$nodes$bfactor)
})

test_that("set_enm stops on a network with negative springs", {
  # kij_hnm is negative below 2390/860 = 2.78 A; with side-chain nodes 2acy has
  # two pairs closer than that, which give the Hessian negative eigenvalues
  prot_ca <- set_enm(pdb_2acy_A, node = "ca", model = "hnm", d_max = 10.5)
  expect_equal(get_nmodes(prot_ca), 3 * get_nsites(prot_ca) - 6)
  expect_error(set_enm(pdb_2acy_A, node = "sc", model = "hnm", d_max = 10.5),
               "negative eigenvalue")
})
