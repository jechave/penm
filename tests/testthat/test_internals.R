# Two kinds of test here, deliberately:
#
#   - CLOSED FORM: the expected value is written down from the definition, not
#     produced by running the function. kij_pfanm(2) is 0.25 because 1/2^2 is
#     0.25. These catch an implementation that stops matching its definition.
#
#   - FROZEN: the expected value comes from fixtures/internals_expected.rda.
#     These catch drift -- any change in any number, including in constants
#     (Hinsen's 860/2390/1280000, REACH's 712/6.92/32.0) that cannot be
#     re-derived here without copying them out of the source.
#
# The kij_* functions matter most: they ARE the model definitions, reached only
# by string dispatch from calculate_enm_graph. If one changed, every downstream
# number would change with it, and the delta_* fixtures could not tell -- they
# were generated through these same functions.
#
# Feed every function the input shape the package actually calls it with. An
# earlier version of this file passed the full reduced cmat to dbhat/rwsip/dh/
# logdet, which no code path does: those take the 3x3 single-site block. The
# full matrix is singular by construction (the six rigid-body modes are dropped
# from the spectrum), so its determinant is zero, logdet has no value on it, and
# the frozen numbers were floating-point residue -- they differed by BLAS and
# broke on CI.

load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "internals_expected.rda"))

wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
dij  <- internals_expected$dij
sdij <- internals_expected$sdij


# closed form ------------------------------------------------------------

test_that("kij_anm is a step function of distance, forcing i,i+1 contacts", {
  # kij = k inside d_max OR when |sdij| == 1, else 0. Default d_max = 10, k = 1.
  expect_equal(penm:::kij_anm(c(5, 12), c(4, 4)), c(1, 0))
  expect_equal(penm:::kij_anm(c(12, 12), c(1, -1)), c(1, 1))  # forced regardless of distance
  expect_equal(penm:::kij_anm(c(5, 12), c(4, 4), d_max = 10, k = 3), c(3, 0))
})

test_that("kij_ming_wall gives k inside the cutoff and a*k for i,i+1", {
  # defaults d_max = 10.5, k = 4.5, a = 42, so forced contacts are 42*4.5 = 189
  expect_equal(penm:::kij_ming_wall(c(5, 12), c(4, 4)), c(4.5, 0))
  expect_equal(penm:::kij_ming_wall(c(5, 12), c(1, -1)), c(189, 189))
})

test_that("kij_pfanm is the inverse square of distance", {
  expect_equal(penm:::kij_pfanm(c(2, 4, 10)), c(0.25, 0.0625, 0.01))
})

test_that("kij_hnm0 is a gaussian in distance", {
  # a * exp(-(dij/c)^2), defaults c = 7.5, a = 1
  expect_equal(penm:::kij_hnm0(0), 1)
  expect_equal(penm:::kij_hnm0(7.5), exp(-1))
  expect_equal(penm:::kij_hnm0(3, c = 3, a = 2), 2 * exp(-1))
})

test_that("kij_gnm and kij_pfgnm are aliases", {
  # The anm/gnm distinction is whether 3D space is considered, not how contacts
  # are defined, so the spring constants coincide.
  expect_identical(penm:::kij_gnm, penm:::kij_anm)
  expect_identical(penm:::kij_pfgnm, penm:::kij_pfanm)
})

test_that("kij_reach returns its fixed constants for the first three separations", {
  expect_equal(penm:::kij_reach(5, sdist = 1), 712)
  expect_equal(penm:::kij_reach(5, sdist = 2), 6.92)
  expect_equal(penm:::kij_reach(5, sdist = 3), 32.0)
})

test_that("distribution measures are degenerate on identical inputs", {
  # 3x3 single-site block, the shape these are called with in the package.
  # Not the full cmat: it is singular, so dbhat's logdet calls are undefined
  # there and rwsip's eigendecomposition produces NaN.
  b <- site_block(get_cmat(wt), 11L)
  expect_equal(penm:::rwsip(b, b), 1)   # perfect similarity
  expect_equal(penm:::dbhat(b, b), 0)   # zero distance
  expect_equal(penm:::dh(b, b), 0)      # zero entropy difference
})

test_that("small utilities match their definitions", {
  m <- matrix(c(1, 7, 7, 2), 2, 2)
  expect_equal(penm:::tr(m), 3)
  expect_equal(penm:::logdet(diag(c(2, 4))), log(8))
  expect_equal(penm:::xyz_indices_site(1), c(1, 2, 3))
  expect_equal(penm:::xyz_indices_site(3), c(7, 8, 9))
  expect_equal(penm:::xyz_indices_site(c(1, 2)), c(1, 2, 3, 4, 5, 6))
  expect_equal(penm:::my_quad_form(c(1, 0), diag(2), c(1, 0)), 1)
  # reduce_matrix sums the trace of each 3x3 block
  expect_equal(penm:::reduce_matrix(diag(6)), matrix(c(3, 0, 0, 3), 2, 2))
})

test_that("my_as_xyz reshapes to 3 rows and rejects a bad length", {
  expect_equal(dim(penm:::my_as_xyz(1:6)), c(3L, 2L))
  expect_error(penm:::my_as_xyz(1:5), "Stopping execution")
})

test_that("v_dij is the spring energy 1/2 k (dij - lij)^2 plus v0ij", {
  expect_equal(penm:::v_dij(dij = 6, v0ij = 0, kij = 2, lij = 4), 4)
  expect_equal(penm:::v_dij(dij = 4, v0ij = 5, kij = 2, lij = 4), 5)
})


# frozen ------------------------------------------------------------------

test_that("kij_* match frozen values", {
  expect_equal(penm:::kij_anm(dij, sdij), internals_expected$kij_anm)
  expect_equal(penm:::kij_gnm(dij, sdij), internals_expected$kij_gnm)
  expect_equal(penm:::kij_hnm(dij), internals_expected$kij_hnm)
  expect_equal(penm:::kij_hnm0(dij), internals_expected$kij_hnm0)
  expect_equal(penm:::kij_ming_wall(dij, sdij), internals_expected$kij_ming_wall)
  expect_equal(penm:::kij_pfanm(dij), internals_expected$kij_pfanm)
  expect_equal(penm:::kij_pfgnm(dij), internals_expected$kij_pfgnm)
})

test_that("kij_reach matches frozen values on every branch", {
  expect_equal(penm:::kij_reach(dij, sdist = 1), internals_expected$kij_reach_sd1)
  expect_equal(penm:::kij_reach(dij, sdist = 2), internals_expected$kij_reach_sd2)
  expect_equal(penm:::kij_reach(dij, sdist = 3), internals_expected$kij_reach_sd3)
  expect_equal(penm:::kij_reach(dij, sdist = 5, same_chain = TRUE), internals_expected$kij_reach_in)
  expect_equal(penm:::kij_reach(dij, sdist = 5, same_chain = FALSE), internals_expected$kij_reach_ex)
})

test_that("distribution measures match frozen values", {
  # Compared against a DIFFERENT ENM model, not against the lfenm mutant: lfenm
  # leaves kmat and hence cmat unchanged, so that comparison would freeze the
  # degenerate 0/1/0 answers and could not detect a broken function.
  wt_anm <- set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = 10.5, frustrated = FALSE)
  bwt  <- site_block(get_cmat(wt),     internals_expected$site_block)
  bmut <- site_block(get_cmat(wt_anm), internals_expected$site_block)
  expect_equal(penm:::dbhat(bwt, bmut), internals_expected$dbhat)
  expect_equal(penm:::rwsip(bwt, bmut), internals_expected$rwsip)
  expect_equal(penm:::dh(bwt, bmut), internals_expected$dh)
  expect_equal(penm:::dbhat(bwt, bmut, normalize = TRUE), internals_expected$dbhat_norm)
  expect_equal(penm:::rwsip(bwt, bmut, normalize = TRUE), internals_expected$rwsip_norm)
  expect_equal(penm:::dh(bwt, bmut, normalize = TRUE), internals_expected$dh_norm)
})

test_that("utilities and energies match frozen values", {
  cwt <- get_reduced_cmat(wt)
  expect_equal(penm:::tr(cwt), internals_expected$tr)   # a plain diagonal sum, fine here
  # logdet on the 3x3 block, not on cwt: cwt is singular (rigid-body modes
  # dropped), so its determinant is zero and log(det) has no value.
  expect_equal(penm:::logdet(site_block(get_cmat(wt), internals_expected$site_block)),
               internals_expected$logdet)
  expect_equal(enm_g_entropy(wt, beta_boltzmann()), internals_expected$enm_g_entropy)
  expect_equal(penm:::v_dij(dij = c(3, 5, 7), v0ij = 0, kij = 2, lij = 4),
               internals_expected$v_dij)
  expect_equal(penm:::cn_xyz(get_xyz(wt), 10.5), internals_expected$cn_xyz)
  expect_equal(penm:::wcn_xyz(get_xyz(wt)), internals_expected$wcn_xyz)
  expect_equal(penm:::cn_graph(get_graph(wt)), internals_expected$cn_graph)
})
