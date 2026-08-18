# Refresh data to use in test_internals.R
#
# Frozen output of the internal helpers: the kij_* spring-constant functions
# that define the ENM model variants, the distribution-comparison measures in
# utils.R, and the small pure utilities.
#
# Same standing as the other fixtures: a REGRESSION anchor, not a correctness
# claim. A mismatch means a number moved -- find out why before regenerating.
#
# These matter more than their size suggests. The kij_* functions ARE the model
# definitions, and they are reached only by string dispatch
# (match.fun(paste0("kij_", model)) in calculate_enm_graph). If one of them
# changed, every downstream number would change together, and the delta_*
# fixtures could not reveal it -- they were generated through the same
# functions.

## ----------------------------------------------------------------------------------------------------------------------
# load libraries
library(tidyverse)
library(bio3d)
library(penm)
library(jefuns)
library(here)


## ----------------------------------------------------------------------------------------------------------------------
# fixed inputs: a spread of distances either side of the usual cutoffs, and
# sequence separations covering the i,i+1 special case.

dij  <- c(0.5, 1, 2, 3.5, 3.999, 4, 4.001, 6, 7.5, 9, 10, 10.5, 11, 15)
sdij <- c(1, -1, 2, 3, 1, 5, -1, 2, 4, 1, 6, -1, 3, 8)

internals_expected <- list(
  dij  = dij,
  sdij = sdij,
  # kij_* : the model definitions. kij_gnm and kij_pfgnm are aliases of
  # kij_anm and kij_pfanm (same contact graph; the anm/gnm distinction is
  # whether 3D space is considered, not how contacts are defined).
  kij_anm       = penm:::kij_anm(dij, sdij),
  kij_gnm       = penm:::kij_gnm(dij, sdij),
  kij_hnm       = penm:::kij_hnm(dij),
  kij_hnm0      = penm:::kij_hnm0(dij),
  kij_ming_wall = penm:::kij_ming_wall(dij, sdij),
  kij_pfanm     = penm:::kij_pfanm(dij),
  kij_pfgnm     = penm:::kij_pfgnm(dij),
  # kij_reach dispatches on a scalar sdist, so record one value per branch
  kij_reach_sd1 = penm:::kij_reach(dij, sdist = 1),
  kij_reach_sd2 = penm:::kij_reach(dij, sdist = 2),
  kij_reach_sd3 = penm:::kij_reach(dij, sdist = 3),
  kij_reach_in  = penm:::kij_reach(dij, sdist = 5, same_chain = TRUE),
  kij_reach_ex  = penm:::kij_reach(dij, sdist = 5, same_chain = FALSE)
)

## ----------------------------------------------------------------------------------------------------------------------
# distribution-comparison measures, on two real covariance matrices

load(here("tests/testthat/fixtures/pdb_2acy_A.rda"))

wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
mut <- get_mutant_site(wt, site_mut = 80, mutation = 1,
                       mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1,
                       seed = 241956)

cwt  <- get_reduced_cmat(wt)
cmut <- get_reduced_cmat(mut)

internals_expected$dbhat        <- penm:::dbhat(cwt, cmut)
internals_expected$rwsip        <- penm:::rwsip(cwt, cmut)
internals_expected$dh           <- penm:::dh(cwt, cmut)
internals_expected$dbhat_norm   <- penm:::dbhat(cwt, cmut, normalize = TRUE)
internals_expected$rwsip_norm   <- penm:::rwsip(cwt, cmut, normalize = TRUE)
internals_expected$dh_norm      <- penm:::dh(cwt, cmut, normalize = TRUE)

# small pure utilities
internals_expected$tr           <- penm:::tr(cwt)
internals_expected$logdet       <- penm:::logdet(cwt)
internals_expected$enm_g_entropy <- enm_g_entropy(wt, beta_boltzmann())
internals_expected$v_dij        <- penm:::v_dij(dij = c(3, 5, 7), v0ij = 0, kij = 2, lij = 4)
internals_expected$cn_xyz       <- penm:::cn_xyz(get_xyz(wt), 10.5)
internals_expected$wcn_xyz      <- penm:::wcn_xyz(get_xyz(wt))
internals_expected$cn_graph     <- penm:::cn_graph(get_graph(wt))

save(internals_expected, file = here("tests/testthat/fixtures/internals_expected.rda"))
