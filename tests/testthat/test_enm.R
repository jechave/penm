library(here)
library(tidyverse)
library(bio3d)
library(jefuns)
library(penm)

load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "prot_2acy_A_ming_wall_ca.rda"))


test_that("set_enm gets prot  equal to prot_2acy_A", {
  expect_equal(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE),
               prot_2acy_A_ming_wall_ca)
})

test_that("set_enm works with beta carbon nodes", {
  # Test that CB nodes can be created
  prot_cb <- set_enm(pdb_2acy_A, node = "cb", model = "anm", d_max = 12.0, frustrated = FALSE)

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
  prot_cb1 <- set_enm(pdb_2acy_A, node = "cb", model = "anm", d_max = 12.0, frustrated = FALSE)
  prot_cb2 <- set_enm(pdb_2acy_A, node = "beta", model = "anm", d_max = 12.0, frustrated = FALSE)

  # Both should produce the same result
  expect_equal(prot_cb1$nodes$xyz, prot_cb2$nodes$xyz)
  expect_equal(prot_cb1$nodes$bfactor, prot_cb2$nodes$bfactor)
})
