library(here)
library(tidyverse)
library(jefuns)
library(bio3d)

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
                       seed = 241956)

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

test_that("delta_motion by site match frozen values", {
  expect_equal(delta_motion_dmsfi(wt, mut), delta_expected$dmsfi)
  expect_equal(delta_motion_dbhati(wt, mut), delta_expected$dbhati)
  expect_equal(delta_motion_rwsipi(wt, mut), delta_expected$rwsipi)
  expect_equal(delta_motion_dhi(wt, mut), delta_expected$dhi)
})

test_that("delta_motion by mode match frozen values", {
  expect_equal(delta_motion_dmsfn(wt, mut), delta_expected$dmsfn)
  expect_equal(delta_motion_dhn(wt, mut), delta_expected$dhn)
  expect_equal(delta_motion_rwsipn(wt, mut), delta_expected$rwsipn)
  expect_equal(delta_motion_nhn(wt, mut), delta_expected$nhn)
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
                          seed = 241956)
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
