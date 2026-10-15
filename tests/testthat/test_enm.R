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

test_that("set_enm allows hnm only with alpha-carbon nodes", {
  # kij_hnm is negative below 2390/860 = 2.78 A; with side-chain or beta-carbon
  # nodes 2acy has pairs closer than that
  prot_ca <- set_enm(pdb_2acy_A, node = "ca", model = "hnm", d_max = 10.5)
  expect_equal(get_nmodes(prot_ca), 3 * get_nsites(prot_ca) - 6)
  expect_error(set_enm(pdb_2acy_A, node = "sc", model = "hnm", d_max = 10.5),
               "C-alpha nodes only")
  expect_error(set_enm(pdb_2acy_A, node = "cb", model = "hnm", d_max = 10.5),
               "C-alpha nodes only")
})

test_that("the normal modes stop on a network with negative springs", {
  # built past set_enm's node check: hnm springs on side-chain nodes, two of
  # them negative, give the Hessian negative eigenvalues
  nodes <- penm:::calculate_enm_nodes(pdb_2acy_A, "sc")
  graph <- penm:::calculate_enm_graph(nodes$xyz, nodes$pdb_site, "hnm", 10.5, list())
  expect_equal(sum(graph$kij < 0), 2)
  eij <- penm:::calculate_enm_eij(nodes$xyz, graph$i, graph$j)
  kmat <- penm:::calculate_enm_kmat(graph, eij, nodes$nsites)
  expect_error(penm:::calculate_enm_nma(kmat), "negative eigenvalue")
})

test_that("fluctuations are covariances at inverse temperature beta, in A^2", {
  wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5)
  # along each mode, kT / lambda_n, with kT in kcal/mol at the default 298 K
  expect_equal(get_msf_mode(wt), 1.987204e-3 * 298 / get_evalue(wt))
  # proportional to temperature
  ratio <- get_msf_site(wt, beta_boltzmann(T = 310)) / get_msf_site(wt, beta_boltzmann(T = 298))
  expect_equal(unname(ratio), rep(310 / 298, get_nsites(wt)))
  # cmat is K+ / beta, with K+ the pseudo-inverse of kmat
  k <- get_kmat(wt)
  expect_equal(k %*% get_cmat(wt) %*% k, k / beta_boltzmann())
})

test_that("a prot from penm < 0.6.0, which stored nma$cmat, is rejected", {
  wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5)
  old <- wt
  names(old$nma)[names(old$nma) == "kmat_pinv"] <- "cmat"
  expect_error(get_cmat(old), "rebuild it with set_enm")
})

test_that("set_enm stores how the protein mutates, with defaults", {
  wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5)
  expect_identical(get_enm_param(wt), list(
    node = "ca", model = "ming_wall", d_max = 10.5, d_max_graph = 10.5, kij_par = list(),
    mut_model = "lfenm", ensemble = 1L, alphabet_size = 10L, mut_dl_sigma = 0.3, mut_sd_min = 2L
  ))

  # the mutation parameters reach the mutant: mut_sd_min = 1 perturbs the
  # i,i+1 bonds as well, which the default leaves alone
  wt_1 <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, mut_sd_min = 1)
  bonds <- get_graph(wt)$sdij == 1
  bonds_changed <- function(p) sum(get_mutant_site(p, 80, 1)$graph$lij[bonds] != p$graph$lij[bonds])
  expect_equal(bonds_changed(wt), 0)
  expect_equal(bonds_changed(wt_1), 2)

  # and so does mut_dl_sigma: the same draws, scaled
  wt_s <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, mut_dl_sigma = 0.6)
  dlij <- function(p) get_mutant_site(p, 80, 1)$graph$lij - p$graph$lij
  expect_equal(dlij(wt_s), 2 * dlij(wt))
})

test_that("set_enm passes further parameters to the spring-constant function", {
  # lfenm with a smoothed cutoff: kij is k(dij) of the smooth function
  wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall_smooth", d_max = 10.5, d_max_width = 1)
  graph <- get_graph(wt)
  expect_identical(get_enm_param(wt)$kij_par, list(d_max_width = 1))
  expect_equal(graph$kij, penm:::kij_ming_wall_smooth(graph$dij, graph$sdij, d_max = 10.5, d_max_width = 1))
  expect_true(any(graph$kij != penm:::kij_ming_wall(graph$dij, graph$sdij, d_max = 10.5)))
})

test_that("set_enm validates its arguments", {
  b <- function(...) set_enm(pdb_2acy_A, node = "ca", d_max = 10.5, ...)
  expect_error(b(model = "ming_wall", mut_model = "nope"), "mut_model must be")
  expect_error(b(model = "ming_wall", d_max_graph = 14), "must equal d_max")
  expect_error(b(model = "ming_wall", d_max_graph = 9), "must not be smaller")
  expect_error(b(model = "ming_wall", wdth = 1), "does not take parameter")
  expect_error(b(model = "ming_wall", 1), "must be named")
  expect_error(b(model = "ming_wall_smooth"), "\"d_max_width\" is missing")
  expect_error(b(model = "no_such_model"), "kij_no_such_model")
  expect_error(b(model = "ming_wall", ensemble = 1.5), "ensemble")
  expect_error(b(model = "ming_wall", alphabet_size = 1), "alphabet_size")
  expect_error(b(model = "ming_wall", mut_dl_sigma = 0), "mut_dl_sigma")
  expect_error(b(model = "ming_wall", mut_sd_min = 0), "mut_sd_min")
})
