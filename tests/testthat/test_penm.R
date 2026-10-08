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


# alleles -------------------------------------------------------------------

wt_alleles <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, ensemble = 3L)

test_that("get_mutant_site validates the site and the allele", {
  expect_error(get_mutant_site(wt_alleles, 0, 1), "site_mut must be")
  expect_error(get_mutant_site(wt_alleles, get_nsites(wt_alleles) + 1, 1), "site_mut must be")
  expect_error(get_mutant_site(wt_alleles, 80, 10), "mutation must be an allele")
  expect_error(get_mutant_site(wt_alleles, 80, -1), "mutation must be an allele")
  expect_error(get_mutant_site(wt_alleles, 80, 1.5), "mutation must be an allele")
  no_sequence <- wt_alleles
  no_sequence$nodes$sequence <- NULL
  expect_error(get_mutant_site(no_sequence, 80, 1), "Rebuild it with set_enm")
})

test_that("the allele a site already has returns the protein unchanged", {
  expect_identical(get_mutant_site(wt_alleles, 80, 0), wt_alleles)
  m <- get_mutant_site(wt_alleles, 80, 3)
  expect_identical(get_mutant_site(m, 80, 3), m)
})

test_that("mutating back restores the wild type", {
  m <- get_mutant_site(wt_alleles, 80, 3)
  back <- get_mutant_site(m, 80, 0)
  expect_gt(max(abs(get_xyz(m) - get_xyz(wt_alleles))), 0.01)
  expect_identical(back$graph$lij, wt_alleles$graph$lij)
  expect_identical(back$nodes$sequence, wt_alleles$nodes$sequence)
  # as.vector: an lfenm mutant's xyz is a 3N x 1 matrix, the wild type's a vector
  expect_equal(as.vector(get_xyz(back)), get_xyz(wt_alleles))
})

test_that("an lfenm mutant is the wild type plus the response to lij - l0ij, whatever the path", {
  # two paths to one sequence, through different intermediate alleles
  path_1 <- get_mutant_site(get_mutant_site(get_mutant_site(wt_alleles, 80, 3), 20, 5), 80, 7)
  path_2 <- get_mutant_site(get_mutant_site(wt_alleles, 80, 7), 20, 5)
  expect_identical(path_1$nodes$sequence, path_2$nodes$sequence)
  expect_identical(path_1$graph$lij, path_2$graph$lij)

  # closed form, in one step from the wild type: kmat and the edge directions
  # are the wild type's, so the response is linear in lij - l0ij
  delta_lij <- path_1$graph$lij - wt_alleles$graph$l0ij
  f <- penm:::calculate_force(wt_alleles, delta_lij)
  xyz_direct <- get_xyz(wt_alleles) + as.vector(get_cmat(wt_alleles) %*% f)
  expect_gt(max(abs(xyz_direct - get_xyz(wt_alleles))), 0.01)
  # as.vector: an lfenm mutant's xyz is a 3N x 1 matrix
  expect_equal(as.vector(get_xyz(path_1)), xyz_direct)
  expect_equal(as.vector(get_xyz(path_2)), xyz_direct)
})

test_that("an allele makes the same change to an edge in lfenm and in genm", {
  # same ensemble, mut_dl_sigma and mut_sd_min; genm has more edges (d_max_graph)
  wt_genm <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, ensemble = 3L,
                     mut_model = "genm", d_max_graph = 14)
  m_lfenm <- get_mutant_site(wt_alleles, 80, 3)
  m_genm <- penm:::genm_mutate(wt_genm, 80, 3)
  change <- function(p) tibble::tibble(edge = p$graph$edge, change = p$graph$lij - p$graph$l0ij)
  shared <- dplyr::inner_join(change(m_lfenm), change(m_genm), by = "edge", suffix = c("_lfenm", "_genm"))
  expect_equal(nrow(shared), nrow(m_lfenm$graph))
  expect_gt(sum(shared$change_lfenm != 0), 10)
  expect_identical(shared$change_lfenm, shared$change_genm)
})
