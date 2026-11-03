# Run from the package root: Rscript dev/explorations/motion_comparisons/motion_identities.R
suppressMessages(devtools::load_all(quiet = TRUE))
wt  <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, mut_model = "genm", d_max_graph = 14)
mut <- set_enm_nma(get_mutant_site(wt, site_mut = 11, mutation = 1))
s2  <- get_msf_mode(wt)
dmsfn <- delta_motion_dmsfn(wt, mut)
# 1. dhn is a transform of dmsfn
cat("dhn == 1/2 log(1 + dmsfn/s2):", isTRUE(all.equal(delta_motion_dhn(wt, mut), 0.5 * log(1 + dmsfn / s2))), "\n")
# 2. rwsipn is a transform of dmsfn, given the paired eigenvalue normaliser
norm <- sum(s2 * get_msf_mode(mut))
cat("rwsipn == sqrt(s2 (s2 + dmsfn) / sum(s2 mu)):", isTRUE(all.equal(delta_motion_rwsipn(wt, mut), sqrt(s2 * (s2 + dmsfn) / norm))), "\n")
# 3. sum of dhn vs the entropy change over own modes
dH_total <- 0.5 * sum(log(get_msf_mode(mut))) - 0.5 * sum(log(s2))
cat("sum(dhn) =", sum(delta_motion_dhn(wt, mut)), " true dH =", dH_total, " -beta*ddg_tds =", -beta_boltzmann() * ddg_tds(wt, mut), "\n")
# 4. KL(mut || wt) on wt's internal subspace, and its per-wt-mode split
U  <- get_umat(wt)
Cb <- t(U) %*% get_cmat(mut) %*% U          # mut covariance in wt mode basis
x  <- diag(Cb) / s2
dr <- as.vector(crossprod(U, get_xyz(mut) - get_xyz(wt)))
kl_mean_n <- 0.5 * dr^2 / s2
kl_diag_n <- 0.5 * (x - 1 - log(x))
mixing    <- 0.5 * (sum(log(diag(Cb))) - as.numeric(determinant(Cb)$modulus))
kl_cov    <- 0.5 * (sum(x) - length(x) + sum(log(s2)) - as.numeric(determinant(Cb)$modulus))
cat("KL cov == sum(diag terms) + mixing:", isTRUE(all.equal(kl_cov, sum(kl_diag_n) + mixing)), "\n")
cat("KL mean part == beta/2 sum(de2n):", isTRUE(all.equal(sum(kl_mean_n), beta_boltzmann() / 2 * sum(delta_structure_de2n(wt, mut)))), "\n")
cat("KL mean =", sum(kl_mean_n), " KL cov diag =", sum(kl_diag_n), " mixing =", mixing, "\n")
cat("range nhn:", range(delta_motion_nhn(wt, mut)), "\n")
cat("smallest eigenvalue gaps (relative):", head(sort(diff(get_evalue(wt)) / head(get_evalue(wt), -1)), 3), "\n")
cat("range rwsipi:", range(delta_motion_rwsipi(wt, mut)), "\n")
