# Regression tests. The expected values are frozen output (see
# fixtures/refresh_test_delta_data.R), not independently verified results: they
# catch a number changing, they do not certify it is right. A failure here means
# something changed -- find out what before touching the fixture.

load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "delta_expected.rda"))
load(test_path("fixtures", "prot_expected.rda"))

wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
              mut_dl_sigma = 0.3, mut_sd_min = 1, ensemble = 1L)
mut <- get_mutant_site(wt, site_mut = 80, mutation = 1)

pdb_site_active <- c(23, 41)


test_that("delta_structure by site match frozen values", {
  kmat_sqrt <- get_kmat_sqrt(wt)
  expect_equal(delta_structure_dr2i(wt, mut), delta_expected$dr2i)
  expect_equal(delta_structure_de2i(wt, mut, kmat_sqrt = kmat_sqrt), delta_expected$de2i)
  expect_equal(delta_structure_df2i(wt, mut), delta_expected$df2i)
})

test_that("delta_structure by mode match frozen values", {
  expect_equal(delta_structure_dr2n(wt, mut), delta_expected$dr2n)
  expect_equal(delta_structure_de2n(wt, mut), delta_expected$de2n)
  expect_equal(delta_structure_df2n(wt, mut), delta_expected$df2n)
})

# The delta_motion_* comparisons cannot be regression-tested against an lfenm
# mutant. lfenm perturbs edge equilibrium lengths and leaves kmat untouched, so
# cmat and the normal modes are IDENTICAL between wt and mutant, and every
# comparison returns its degenerate value whatever the code does. Freezing
# those numbers records values that hold no matter how broken the function is.
# That is how a NaN in the old delta_motion_nhn survived a passing suite: it was
# frozen into the fixture, and expect_equal(NaN, NaN) is TRUE.
#
# So they are tested on a genm mutant, whose network and fluctuations do
# change, against closed forms and against independent computations of the
# same quantity.

wt_genm <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                   mut_model = "genm", d_max_graph = 14)
mut_genm <- set_enm_nma(get_mutant_site(wt_genm, site_mut = 11, mutation = 1))

# mut_genm, rotated by `degrees` about the z axis through the wild type's
# centroid, with its covariance and modes rotated with it
rotate_mutant <- function(degrees) {
  angle <- degrees * pi / 180
  rotation <- matrix(c(cos(angle), sin(angle), 0,
                       -sin(angle), cos(angle), 0,
                       0, 0, 1), nrow = 3)
  wt_xyz <- matrix(get_xyz(wt_genm), nrow = 3)   # column k is site k
  centroid <- rowMeans(wt_xyz)
  wt_xyz_rotated <- rotation %*% (wt_xyz - centroid) + centroid
  superpose_prot(mut_genm, as.vector(wt_xyz_rotated))
}

test_that("lfenm leaves the network, and so the motion, unchanged", {
  expect_equal(get_kmat(wt), get_kmat(mut))
  expect_equal(get_cmat(wt), get_cmat(mut))
  expect_equal(get_evalue(wt), get_evalue(mut))
})

test_that("get_h_site is the entropy of each site's three-dimensional Gaussian", {
  # independently of the determinant: from the variances along the principal axes
  cmat <- get_cmat(wt_genm)
  expected <- sapply(seq(get_nsites(wt_genm)), function(i) {
    variances <- eigen(site_block(cmat, i), symmetric = TRUE)$values
    sum(0.5 * log(2 * pi * exp(1) * variances))
  })
  expect_equal(get_h_site(wt_genm), expected)
})

test_that("get_h_mode along a protein's own modes follows from its eigenvalues", {
  # get_h_mode projects the covariance; get_msf_mode comes from the eigenvalues
  expect_equal(get_h_mode(wt_genm), 0.5 * log(2 * pi * exp(1) * get_msf_mode(wt_genm)))
})

test_that("the size part of dlcov is the entropy change, unsigned and rescaled", {
  # size comes from eigenvalues, dh_site from determinants
  dlcov <- delta_motion_dlcov_site(wt_genm, mut_genm)
  expect_equal(dlcov$size, abs(delta_motion_dh_site(wt_genm, mut_genm)) / sqrt(3))
})

test_that("size, shape and orientation split the log-Euclidean distance exactly", {
  # independently of the split: the Frobenius distance between the matrices
  # 1/2 log C of each site
  half_log <- function(block) {
    axes <- eigen(block, symmetric = TRUE)
    axes$vectors %*% diag(0.5 * log(axes$values)) %*% t(axes$vectors)
  }
  cmat_wt <- get_cmat(wt_genm)
  cmat_mut <- get_cmat(mut_genm)
  expected <- sapply(seq(get_nsites(wt_genm)), function(i) {
    difference <- half_log(site_block(cmat_mut, i)) - half_log(site_block(cmat_wt, i))
    sqrt(sum(difference^2))
  })
  dlcov <- delta_motion_dlcov_site(wt_genm, mut_genm)
  expect_equal(dlcov$dlcov, expected)
  # every part contributes somewhere, so none of them is trivially zero
  expect_gt(max(dlcov$size), 0.01)
  expect_gt(max(dlcov$shape), 0.01)
  expect_gt(max(dlcov$orientation), 0.01)
})

# wt_genm with every site's 3 x 3 covariance set to `block`, and no correlation
# between sites: a protein whose fluctuations are known exactly. get_cmat()
# returns kmat_pinv / beta.
with_site_covariance <- function(block) {
  prot <- wt_genm
  prot$nma$kmat_pinv <- beta_boltzmann() * kronecker(diag(get_nsites(prot)), block)
  prot
}

test_that("dlcov_site gives the closed forms for a change in size, shape or orientation alone", {
  nsites <- get_nsites(wt_genm)
  sphere <- diag(3)
  # standard deviations 4:1:1, with the volume of the sphere
  needle_x <- diag(c(16, 1, 1)) / 16^(1 / 3)
  # the same needle, rotated by 90 degrees about z
  needle_y <- diag(c(1, 16, 1)) / 16^(1 / 3)

  # size alone: every standard deviation doubled, so ln(sigma) grows by ln 2 on each axis
  size_only <- delta_motion_dlcov_site(with_site_covariance(sphere), with_site_covariance(4 * sphere))
  expect_equal(size_only$size, rep(sqrt(3) * log(2), nsites))
  expect_equal(size_only$shape, rep(0, nsites))
  expect_equal(size_only$orientation, rep(0, nsites))
  expect_equal(delta_motion_dh_site(with_site_covariance(sphere), with_site_covariance(4 * sphere)),
               rep(3 * log(2), nsites))

  # shape alone: ln(sigma) minus its mean goes from 0 to ln 4 * (2/3, -1/3, -1/3);
  # a sphere has no orientation
  shape_only <- delta_motion_dlcov_site(with_site_covariance(sphere), with_site_covariance(needle_x))
  expect_equal(shape_only$size, rep(0, nsites))
  expect_equal(shape_only$shape, rep(log(4) * sqrt(2 / 3), nsites))
  expect_equal(shape_only$orientation, rep(0, nsites))

  # orientation alone: elongation ln(sigma_1 / sigma_2) = ln 4 in both, rotated by
  # 90 degrees, so orientation = sqrt(2 * ln(4)^2 * sin(90 degrees)^2)
  orientation_only <- delta_motion_dlcov_site(with_site_covariance(needle_x), with_site_covariance(needle_y))
  expect_equal(orientation_only$size, rep(0, nsites))
  expect_equal(orientation_only$shape, rep(0, nsites))
  expect_equal(orientation_only$orientation, rep(sqrt(2) * log(4), nsites))
})

test_that("a covariance that is not positive definite stops, rather than giving NaN or -Inf", {
  flat <- diag(c(1, 1, 0))
  expect_error(get_h_site(with_site_covariance(flat)), "not positive definite")
  expect_error(delta_motion_dlcov_site(wt_genm, with_site_covariance(flat)), "not positive definite")
})

test_that("a reference mode partly outside the protein's internal space warns, mostly outside stops", {
  # a translation, which no internal mode contains
  translation_x <- rep(c(1, 0, 0), get_nsites(wt_genm))
  translation_x <- translation_x / sqrt(sum(translation_x^2))
  # wt_genm with its first mode tilted towards the translation, so that the
  # reference's first mode keeps only `share` of itself inside the internal space
  with_first_mode_inside <- function(share) {
    prot <- wt_genm
    prot$nma$umat[, 1] <- sqrt(share) * get_umat(wt_genm)[, 1] + sqrt(1 - share) * translation_x
    prot
  }
  expect_warning(delta_motion_nmix_mode(wt_genm, with_first_mode_inside(0.95)), "inside the protein's internal space")
  expect_error(delta_motion_nmix_mode(wt_genm, with_first_mode_inside(0.5)), "too different")
})

test_that("dh_mode takes the variance along the part of each mode inside the mutant's internal space", {
  # independently: project wt's modes onto mut's internal space, explicitly
  wt_modes <- get_umat(wt_genm)
  mut_modes <- get_umat(mut_genm)
  projector <- mut_modes %*% t(mut_modes)
  directions <- projector %*% wt_modes
  lengths <- sqrt(colSums(directions^2))
  directions <- directions %*% diag(1 / lengths)
  variance_mut <- diag(t(directions) %*% get_cmat(mut_genm) %*% directions)
  expected <- 0.5 * log(variance_mut / get_msf_mode(wt_genm))
  expect_equal(delta_motion_dh_mode(wt_genm, mut_genm), expected)
})

test_that("nmix_mode is 1 for a protein against itself, and between 1 and nmodes otherwise", {
  nmodes <- get_nmodes(wt_genm)
  expect_equal(delta_motion_nmix_mode(wt_genm, wt_genm), rep(1, nmodes))
  nmix <- delta_motion_nmix_mode(wt_genm, mut_genm)
  expect_true(all(nmix >= 1 - 1e-12))
  expect_true(all(nmix <= nmodes))
  # the mutant does mix modes, so the bounds above are not met trivially
  expect_gt(max(nmix), 2)
  # exp of the Shannon entropy is never below the inverse participation ratio,
  # 1 / sum(p^2), of the same weights (the Renyi entropies decrease with order)
  overlap_squared <- crossprod(get_umat(wt_genm), get_umat(mut_genm))^2
  weight <- overlap_squared / rowSums(overlap_squared)
  inverse_participation <- 1 / rowSums(weight^2)
  expect_true(all(nmix >= inverse_participation * (1 - 1e-12)))
})

test_that("a protein compared with itself gives zero change and zero dlcov", {
  nsites <- get_nsites(wt_genm)
  nmodes <- get_nmodes(wt_genm)
  expect_equal(delta_motion_dh_site(wt_genm, wt_genm), rep(0, nsites))
  self <- delta_motion_dlcov_site(wt_genm, wt_genm)
  expect_equal(self$size, rep(0, nsites))
  expect_equal(self$shape, rep(0, nsites))
  # rounding in the axis overlaps (~1e-16) comes out of the square root as ~1e-8
  expect_lt(max(self$orientation), 1e-7)
  expect_equal(delta_motion_dh_mode(wt_genm, wt_genm), rep(0, nmodes))
})

test_that("comparisons that use directions stop on proteins that are not superposed", {
  rotated <- rotate_mutant(1)
  expect_error(delta_motion_dlcov_site(wt_genm, rotated), "not superposed")
  expect_error(delta_motion_dh_mode(wt_genm, rotated), "not superposed")
  expect_error(delta_motion_nmix_mode(wt_genm, rotated), "not superposed")
  expect_error(get_h_mode(rotated, reference = wt_genm), "not superposed")
  # a slight rotation warns instead
  expect_warning(delta_motion_dlcov_site(wt_genm, rotate_mutant(0.01)), "not exactly superposed")
  # the volume of a site's fluctuations does not depend on the frame
  expect_equal(delta_motion_dh_site(wt_genm, rotated), delta_motion_dh_site(wt_genm, mut_genm))
})

test_that("motion comparisons return finite values, never NaN", {
  dlcov <- delta_motion_dlcov_site(wt_genm, mut_genm)
  values <- list(
    dh_site     = delta_motion_dh_site(wt_genm, mut_genm),
    dlcov       = dlcov$dlcov,
    size        = dlcov$size,
    shape       = dlcov$shape,
    orientation = dlcov$orientation,
    dh_mode     = delta_motion_dh_mode(wt_genm, mut_genm),
    nmix_mode   = delta_motion_nmix_mode(wt_genm, mut_genm)
  )
  for (name in names(values)) {
    expect_true(all(is.finite(values[[name]])), label = name)
  }
})

test_that("delta_energy functions match frozen values", {
  expect_equal(ddg_dv(wt, mut), delta_expected$ddg_dv)
  expect_equal(ddg_tds(wt, mut), delta_expected$ddg_tds)
  expect_equal(delta_energy_dvs(wt, mut), delta_expected$delta_energy_dvs)
  expect_equal(ddgact_dv(wt, mut, pdb_site_active = pdb_site_active), delta_expected$ddgact_dv)
  expect_equal(ddgact_tds(wt, mut, pdb_site_active = pdb_site_active), delta_expected$ddgact_tds)
})

test_that("prot profiles match frozen values", {
  expect_equal(get_cn(wt), prot_expected$cn)
  expect_equal(get_wcn(wt), prot_expected$wcn)
  expect_equal(get_msf_site(wt), prot_expected$msf_site)
  expect_equal(get_mlms(wt), prot_expected$mlms)
  expect_equal(get_msf_mode(wt), prot_expected$msf_mode)
  expect_equal(get_dactive(wt, pdb_site_active), prot_expected$dactive)
})

test_that("prot matrices match frozen values", {
  expect_equal(get_umat2(wt), prot_expected$umat2)
  expect_equal(get_msf_site_mode(wt), prot_expected$msf_site_mode)
  expect_equal(get_rho_matrix(wt), prot_expected$rho_matrix)
  expect_equal(get_reduced_cmat(wt), prot_expected$reduced_cmat)
  expect_equal(get_reduced_kmat(wt), prot_expected$reduced_kmat)
  expect_equal(get_kmat_sqrt(wt), prot_expected$kmat_sqrt)
  expect_equal(get_cmat_sqrt(wt), prot_expected$cmat_sqrt)
})

test_that("prot accessors and energies match frozen values", {
  expect_equal(get_nsites(wt), prot_expected$nsites)
  expect_equal(get_site(wt), prot_expected$site)
  expect_equal(get_pdb_site(wt), prot_expected$pdb_site)
  expect_equal(get_xyz(wt), prot_expected$xyz)
  expect_equal(get_bfactor(wt), prot_expected$bfactor)
  expect_equal(get_enm_param(wt), prot_expected$enm_param)
  expect_equal(dgact_dv(wt, ideal = wt, pdb_site_active = pdb_site_active), prot_expected$dgact_dv)
  expect_equal(dgact_tds(wt, ideal = wt, pdb_site_active = pdb_site_active), prot_expected$dgact_tds)
  expect_equal(enm_v_min(wt), prot_expected$enm_v_min)
})

test_that("dr2i and dr2n are the same displacement in two bases", {
  # dr2i is per site, dr2n per mode; they decompose one displacement, so the
  # totals must agree. Same for the deformation energy de2i / de2n.
  expect_equal(sum(delta_structure_dr2i(wt, mut)), sum(delta_structure_dr2n(wt, mut)))
  expect_equal(sum(delta_structure_de2i(wt, mut, kmat_sqrt = get_kmat_sqrt(wt))),
               sum(delta_structure_de2n(wt, mut)))
})

test_that("mutation = 0 leaves every response at zero", {
  # get_mutant_site() returns wt unmutated when mutation is 0, so every
  # wt-vs-mutant difference must vanish.
  mut0 <- get_mutant_site(wt, site_mut = 80, mutation = 0)
  nsites <- get_nsites(wt)
  nmodes <- length(get_msf_mode(wt))
  expect_equal(delta_structure_dr2i(wt, mut0), rep(0, nsites))
  expect_equal(delta_structure_df2i(wt, mut0), rep(0, nsites))
  expect_equal(delta_structure_dr2n(wt, mut0), rep(0, nmodes))
  expect_equal(delta_structure_df2n(wt, mut0), rep(0, nmodes))
  expect_equal(delta_motion_dh_site(wt, mut0), rep(0, nsites))
  expect_equal(delta_motion_dh_mode(wt, mut0), rep(0, nmodes))
  expect_equal(ddg_dv(wt, mut0), 0)
  expect_equal(delta_energy_dvs(wt, mut0), 0)
})

test_that("the stress energy takes each spring's length in ideal's conformation, by spring", {
  # prot and ideal with different sets of springs: an lfenm protein (springs
  # within d_max = 10.5) against a genm ideal (springs within d_max_graph = 14),
  # so ideal's graph has springs prot lacks, and they interleave with prot's
  ideal <- get_mutant_site(set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
                                   mut_model = "genm", d_max_graph = 14),
                           site_mut = 40, mutation = 2)
  prot <- mut
  expect_gt(nrow(get_graph(ideal)), nrow(get_graph(prot)))
  expect_true(all(get_graph(prot)$edge %in% get_graph(ideal)$edge))

  # independently: each of prot's springs, at its length in ideal's structure
  g <- get_graph(prot)
  d_ideal <- penm:::dij_edge(get_xyz(ideal), g$i, g$j)
  expected <- sum(g$v0ij + 0.5 * g$kij * (d_ideal - g$lij)^2)
  expect_equal(penm:::calculate_vs(prot, ideal), expected)

  # and the other way round: prot has springs ideal lacks, whose lengths come
  # from ideal's coordinates
  g <- get_graph(ideal)
  d_prot <- penm:::dij_edge(get_xyz(prot), g$i, g$j)
  expected <- sum(g$v0ij + 0.5 * g$kij * (d_prot - g$lij)^2)
  expect_equal(penm:::calculate_vs(ideal, prot), expected)
})
