# Regression tests. The expected values are frozen output (see
# fixtures/refresh_test_delta_data.R), not independently verified results: they
# catch a number changing, they do not certify it is right. A failure here means
# something changed -- find out what before touching the fixture.

load(test_path("fixtures", "pdb_2acy_A.rda"))
load(test_path("fixtures", "delta_expected.rda"))
load(test_path("fixtures", "prot_expected.rda"))

wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
mut <- get_mutant_site(wt, site_mut = 80, mutation = 1,
                       mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1,
                       ensemble = 1L)

pdb_site_active <- c(23, 41)


test_that("delta_structure by site match frozen values", {
  kmat_sqrt <- get_kmat_sqrt(wt)
  expect_equal(delta_structure_dr2i(wt, mut), delta_expected$dr2i)
  expect_equal(delta_structure_de2i(wt, mut, kmat_sqrt = kmat_sqrt), delta_expected$de2i)
  expect_equal(delta_structure_df2i(wt, mut), delta_expected$df2i)
  expect_equal(delta_structure_dvmi(wt, mut), delta_expected$dvmi)
  expect_equal(delta_structure_dvsi(wt, mut), delta_expected$dvsi)
  expect_equal(delta_structure_dvsi_same_topology(wt, mut), delta_expected$dvsi_same_topology)
})

test_that("delta_structure by mode match frozen values", {
  expect_equal(delta_structure_dr2n(wt, mut), delta_expected$dr2n)
  expect_equal(delta_structure_de2n(wt, mut), delta_expected$de2n)
  expect_equal(delta_structure_df2n(wt, mut), delta_expected$df2n)
})

# The delta_motion_* measures cannot be regression-tested against an lfenm
# mutant. lfenm perturbs edge equilibrium lengths and leaves kmat untouched, so
# cmat and the normal modes are IDENTICAL between wt and mutant, and every
# motion measure returns its degenerate value. Freezing those numbers records
# 0 / 1 / 0 -- values that hold no matter how broken the function is. That is
# how a NaN in delta_motion_nhn survived a passing suite: it was frozen into
# the fixture, and expect_equal(NaN, NaN) is TRUE.
#
# So assert the invariance itself, in closed form, instead of freezing its
# numeric shadow. These fail on a NaN, which the frozen version could not.
# Exercising the measures on genuinely different networks needs a mutant model
# that rebuilds the contact map -- see dev/ideas.md, to be settled with sclfenm.

test_that("lfenm leaves the network, and so the motion, unchanged", {
  expect_equal(get_kmat(wt), get_kmat(mut))
  expect_equal(get_cmat(wt), get_cmat(mut))
  expect_equal(get_evalue(wt), get_evalue(mut))
})

test_that("delta_motion by site is degenerate under lfenm", {
  nsites <- get_nsites(wt)
  expect_equal(delta_motion_dmsfi(wt, mut), rep(0, nsites))
  expect_equal(delta_motion_dbhati(wt, mut), rep(0, nsites))
  expect_equal(delta_motion_rwsipi(wt, mut), rep(1, nsites))
  expect_equal(delta_motion_dhi(wt, mut), rep(0, nsites))
})

test_that("delta_motion by mode is degenerate under lfenm", {
  nmodes <- get_nmodes(wt)
  # dmsfn and dhn come out as rounding noise around zero rather than exact 0
  expect_equal(delta_motion_dmsfn(wt, mut), rep(0, nmodes), tolerance = 1e-12)
  expect_equal(delta_motion_dhn(wt, mut), rep(0, nmodes), tolerance = 1e-12)
  # rwsipn is weighted by wmat/tr(wmat), a GLOBAL normaliser, so it is not 1
  # per mode even for identical inputs. With overlap = I the weight matrix is
  # diagonal and mode n reduces to msf_n / sqrt(sum(msf^2)) -- so the vector of
  # rwsipn values is a unit vector.
  s2 <- get_msf_mode(wt)
  expect_equal(delta_motion_rwsipn(wt, mut), s2 / sqrt(sum(s2^2)))
  expect_equal(sum(delta_motion_rwsipn(wt, mut)^2), 1)
  # nhn = exp(entropy of the mode-overlap distribution); a perfect one-to-one
  # correspondence is zero entropy, so nhn is 1. This is the assertion the
  # frozen fixture could not make: it caught 29 NaN here.
  expect_equal(delta_motion_nhn(wt, mut), rep(1, nmodes))
})

test_that("motion measures return finite values, never NaN", {
  for (f in c(delta_motion_dmsfi, delta_motion_dbhati, delta_motion_rwsipi,
              delta_motion_dhi, delta_motion_dmsfn, delta_motion_dhn,
              delta_motion_rwsipn, delta_motion_nhn)) {
    v <- f(wt, mut)
    expect_false(any(is.nan(v)))
    expect_true(all(is.finite(v)))
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
  expect_equal(get_stress(wt), prot_expected$stress)
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
  mut0 <- get_mutant_site(wt, site_mut = 80, mutation = 0,
                          mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1,
                          ensemble = 1L)
  nsites <- get_nsites(wt)
  nmodes <- length(get_msf_mode(wt))
  expect_equal(delta_structure_dr2i(wt, mut0), rep(0, nsites))
  expect_equal(delta_structure_df2i(wt, mut0), rep(0, nsites))
  expect_equal(delta_structure_dr2n(wt, mut0), rep(0, nmodes))
  expect_equal(delta_structure_df2n(wt, mut0), rep(0, nmodes))
  expect_equal(delta_motion_dmsfi(wt, mut0), rep(0, nsites))
  expect_equal(delta_motion_dmsfn(wt, mut0), rep(0, nmodes))
  expect_equal(ddg_dv(wt, mut0), 0)
  expect_equal(delta_energy_dvs(wt, mut0), 0)
})
