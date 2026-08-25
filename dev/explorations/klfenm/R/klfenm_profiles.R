## Profiles: how two states differ SITE BY SITE and MODE BY MODE.
##
## Scalars hide what matters here. A correlation of 0.98 between two RMSF
## profiles is compatible with one site being wrong by 25%; an RMSIP of 0.95
## over a 10-mode subspace is compatible with individual modes being unrelated.
## So everything below returns a vector, and the scalar summaries are computed
## from it rather than instead of it.

## ------------------------------------------------------- structure by site --

## Squared displacement of each site between two structures, after removing the
## rigid-body part. Without the superposition this measures how far the molecule
## drifted, not how it changed shape.
superpose <- function(X, Y) {
  cX <- rowMeans(X); cY <- rowMeans(Y)
  A <- X - cX; B <- Y - cY
  s <- svd(A %*% t(B))
  dsign <- sign(det(s$v %*% t(s$u)))
  t(s$v %*% diag(c(1, 1, dsign)) %*% t(s$u)) %*% B
}

klfenm_dr2_site <- function(ref, mut) {
  B <- superpose(ref$R, mut$R)
  colSums((B - (ref$R - rowMeans(ref$R)))^2)
}

klfenm_rmsd <- function(ref, mut) sqrt(mean(klfenm_dr2_site(ref, mut)))

## ------------------------------------------------------------ motion by site --

## RMSF per site, sqrt of the mean-square fluctuation.
klfenm_rmsf_site <- function(st, beta = 1, frustrated = TRUE) {
  sqrt(klfenm_msf_site(klfenm_cmat(klfenm_kmat(st, frustrated = frustrated), beta)))
}

## ------------------------------------------------------------ modes by mode --

## Compare two spectra mode by mode.
##
## Modes are matched by MAXIMUM OVERLAP, not by index: when two eigenvalues are
## close the eigenvectors can swap order, and index-matching then compares
## unrelated modes and reports a spuriously low overlap. `n_reordered` counts
## how often the two disagree, so the reader can see whether it mattered.
klfenm_mode_comparison <- function(nma_a, nma_b, nmodes = 20) {
  n <- min(nmodes, ncol(nma_a$vector), ncol(nma_b$vector))
  Ua <- nma_a$vector[, seq_len(n), drop = FALSE]
  Ub <- nma_b$vector[, seq_len(n), drop = FALSE]
  ov <- abs(crossprod(Ua, Ub))                    # |<u_a, u_b>|

  best  <- apply(ov, 1, which.max)                # greedy match by overlap
  ov_by_index   <- diag(ov)
  ov_by_overlap <- ov[cbind(seq_len(n), best)]

  la <- nma_a$value[seq_len(n)]; lb <- nma_b$value[seq_len(n)]
  list(
    n_modes        = n,
    eigen_a        = la,
    eigen_b        = lb,
    eigen_rel_diff = (lb - la) / la,
    overlap_index  = ov_by_index,
    overlap_best   = ov_by_overlap,
    matched_to     = best,
    n_reordered    = sum(best != seq_len(n)),
    rmsip          = sqrt(sum(ov^2) / n),
    overlap_matrix = ov
  )
}

## RWSIP-style subspace similarity, weighted by 1/eigenvalue so the soft modes
## (which dominate the fluctuations) count most.
klfenm_rwsip <- function(nma_a, nma_b, nmodes = 20) {
  n <- min(nmodes, ncol(nma_a$vector), ncol(nma_b$vector))
  ov <- crossprod(nma_a$vector[, seq_len(n), drop = FALSE],
                  nma_b$vector[, seq_len(n), drop = FALSE])
  w  <- 1 / nma_a$value[seq_len(n)]
  sqrt(sum(w * rowSums(ov^2)) / sum(w))
}

## ------------------------------------------------------------- convenience --

## Everything about one state that a comparison needs, computed once.
klfenm_observe <- function(st, beta = 1, frustrated = TRUE, nmodes = 20) {
  K <- klfenm_kmat(st, frustrated = frustrated)
  s <- klfenm_nma(K)
  list(state = st, kmat = K, nma = s,
       rmsf  = sqrt(klfenm_msf_site(klfenm_cmat(K, beta))),
       ts    = klfenm_entropy(K, beta),
       n_zero = s$n_zero,
       lowest_raw = min(s$raw_value),
       energy = klfenm_energy(st),
       strain = klfenm_strain(st),
       n_active = n_active(st))
}

## Relative difference profile, in %, with the worst offenders named.
## Reports the whole vector; `worst` is a reading aid, not a substitute.
rel_profile <- function(a, b, top = 5) {
  rel <- 100 * (b - a) / a
  o <- order(abs(rel), decreasing = TRUE)[seq_len(min(top, length(rel)))]
  list(rel = rel, mean = mean(rel), sd = sd(rel),
       min = min(rel), max = max(rel),
       cor = stats::cor(a, b),
       worst_idx = o, worst_val = rel[o])
}
