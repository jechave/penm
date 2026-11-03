# Run from the package root: Rscript dev/explorations/motion_comparisons/metrics_scan.R
#
# Distribution metrics between the wild type and each single-point genm mutant,
# each decomposed into per-site and per-mode contributions.
#
# Both distributions are Gaussians on the wild type's internal subspace (the mutant
# is superposed onto the wild type by genm). Work in wt mode coordinates:
#   a  = wt variances along wt modes, sigma_n^2 = 1 / (beta lambda_n)  (diagonal)
#   b  = mutant covariance projected onto wt modes, U' C_mut U
#   dr = mean displacement along wt modes, U' (r_mut - r_wt)
# A contribution matrix M in mode coordinates is attributed to mode n by M[n, n]
# and to site i by the trace of site i's 3 x 3 block of U M U'.
out_dir <- "dev/explorations/motion_comparisons"
suppressMessages(devtools::load_all(quiet = TRUE))

wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, mut_model = "genm", d_max_graph = 14,
              ensemble = 1, alphabet_size = 10, mut_dl_sigma = 0.3, mut_sd_min = 2)
u_wt <- get_umat(wt)
beta <- beta_boltzmann()
wt_variance <- 1 / (beta * get_evalue(wt))
xyz_wt <- get_xyz(wt)
nmodes <- length(wt_variance)
nsites <- get_nsites(wt)
alphabet <- 0:(get_enm_param(wt)$alphabet_size - 1)
mutations <- expand.grid(letter = alphabet, site = get_site(wt))
mutations <- mutations[mutations$letter != wt$nodes$sequence[mutations$site], ]

# --- helpers ---------------------------------------------------------------

symmetric_function <- function(m, f) {
  e <- eigen((m + t(m)) / 2, symmetric = TRUE)
  e$vectors %*% (f(e$values) * t(e$vectors))
}
matrix_sqrt <- function(m) symmetric_function(m, function(x) sqrt(pmax(x, 0)))

# per-site sums of a 3N vector
sum_by_site <- function(x3n) colSums(matrix(x3n, nrow = 3))

# attribution of a contribution matrix given in wt mode coordinates
by_mode <- function(m) diag(m)
by_site <- function(m) {
  cartesian_diagonal <- rowSums((u_wt %*% m) * u_wt)
  sum_by_site(cartesian_diagonal)
}
# attribution of a squared vector given in wt mode coordinates
vector_by_mode <- function(z) z^2
vector_by_site <- function(z) sum_by_site(as.vector(u_wt %*% z)^2)

# divergences of the covariance part, as functions of the generalised eigenvalues nu
g_kl <- function(nu) 0.5 * (nu - 1 - log(nu))
g_jeffreys <- function(nu) 0.5 * (nu + 1 / nu - 2)
g_bhattacharyya <- function(nu) 0.5 * log((1 + nu) / (2 * sqrt(nu)))
g_riemann <- function(nu) log(nu)^2

site_block <- function(i) 3 * i - 2:0

# --- one mutant --------------------------------------------------------------

compare_distributions <- function(row) {
  mut <- set_enm_nma(get_mutant_site(wt, mutations$site[row], mutations$letter[row]))
  a <- wt_variance
  b <- crossprod(u_wt, get_cmat(mut, beta) %*% u_wt)
  b <- (b + t(b)) / 2
  dr <- as.vector(crossprod(u_wt, get_xyz(mut) - xyz_wt))
  site <- list()
  mode <- list()

  # Wasserstein-2: mean part |dr|^2; covariance part E|x - T x|^2 under the optimal map T
  site$w2_mean <- vector_by_site(dr)
  mode$w2_mean <- vector_by_mode(dr)
  a_half <- diag(sqrt(a))
  a_inv_half <- diag(1 / sqrt(a))
  middle <- matrix_sqrt(a_half %*% b %*% a_half)
  transport <- a_inv_half %*% middle %*% a_inv_half
  residual <- diag(nmodes) - transport
  w2_cov_matrix <- residual %*% diag(a) %*% t(residual)
  site$w2_cov <- by_site(w2_cov_matrix)
  mode$w2_cov <- by_mode(w2_cov_matrix)
  w2_cov_total <- sum(a) + sum(diag(b)) - 2 * sum(diag(middle))
  stopifnot(isTRUE(all.equal(sum(mode$w2_cov), w2_cov_total)))

  # spectral divergences: generalised eigenvalues of (b, a), attributed with the whitened basis
  whitened <- a_inv_half %*% b %*% a_inv_half
  e <- eigen((whitened + t(whitened)) / 2, symmetric = TRUE)
  nu <- e$values
  stopifnot(all(nu > 0))
  for (name in c("kl", "jeffreys", "bhattacharyya", "riemann")) {
    g <- get(paste0("g_", name))
    contribution <- e$vectors %*% (g(nu) * t(e$vectors))
    site[[paste0(name, "_cov")]] <- by_site(contribution)
    mode[[paste0(name, "_cov")]] <- by_mode(contribution)
  }
  kl_cov_total <- 0.5 * (sum(diag(whitened)) - nmodes - as.numeric(determinant(whitened)$modulus))
  stopifnot(isTRUE(all.equal(sum(mode$kl_cov), kl_cov_total)))

  # KL(mut || wt) mean part: 1/2 |a^(-1/2) dr|^2
  z <- dr / sqrt(a)
  site$kl_mean <- 0.5 * vector_by_site(z)
  mode$kl_mean <- 0.5 * vector_by_mode(z)

  # Frobenius distance of the covariances (covariance part only, A^4)
  delta_c <- b - diag(a)
  mode$frobenius_cov <- rowSums(delta_c^2)
  delta_c_cartesian <- u_wt %*% delta_c %*% t(u_wt)
  site$frobenius_cov <- sum_by_site(rowSums(delta_c_cartesian^2))

  # sums of divergences between marginals: 3 x 3 site blocks, and 1D along wt modes
  c_wt <- u_wt %*% diag(a) %*% t(u_wt)
  c_mut <- u_wt %*% b %*% t(u_wt)
  dr_cartesian <- as.vector(u_wt %*% dr)
  site$marginal_kl <- numeric(nsites)
  site$marginal_w2 <- numeric(nsites)
  for (i in seq_len(nsites)) {
    k <- site_block(i)
    ca <- c_wt[k, k]
    cb <- c_mut[k, k]
    dri <- dr_cartesian[k]
    ca_inv <- solve(ca)
    site$marginal_kl[i] <- 0.5 * (sum(diag(ca_inv %*% cb)) - 3 + log(det(ca) / det(cb))) +
      0.5 * sum(dri * (ca_inv %*% dri))
    ca_half <- matrix_sqrt(ca)
    site$marginal_w2[i] <- sum(dri^2) + sum(diag(ca)) + sum(diag(cb)) -
      2 * sum(diag(matrix_sqrt(ca_half %*% cb %*% ca_half)))
  }
  ratio <- diag(b) / a
  mode$marginal_kl <- g_kl(ratio) + 0.5 * dr^2 / a
  mode$marginal_w2 <- (sqrt(diag(b)) - sqrt(a))^2 + dr^2

  # penm's measures, squared as in the genm vignette
  site$dmsf2 <- delta_motion_dmsfi(wt, mut, beta)^2
  mode$dmsf2 <- delta_motion_dmsfn(wt, mut, beta)^2

  # every exact decomposition must add up to the same total by site and by mode
  exact <- c("w2_mean", "w2_cov", "kl_mean", "kl_cov", "jeffreys_cov", "bhattacharyya_cov",
             "riemann_cov", "frobenius_cov")
  for (name in exact) stopifnot(isTRUE(all.equal(sum(site[[name]]), sum(mode[[name]]))))

  list(site = site, mode = mode, nu_range = range(nu))
}

# --- check one mutant against penm's de2i: KL mean part by site == beta/2 de2i ---
mut_check <- set_enm_nma(get_mutant_site(wt, 11, 1))
dr_check <- as.vector(crossprod(u_wt, get_xyz(mut_check) - xyz_wt))
kl_mean_site <- 0.5 * vector_by_site(dr_check / sqrt(wt_variance))
de2i <- delta_structure_de2i(wt, mut_check, kmat_sqrt = get_kmat_sqrt(wt))
cat("KL mean part by site == beta/2 de2i:", isTRUE(all.equal(kl_mean_site, beta / 2 * de2i)), "\n")

# --- the scan ------------------------------------------------------------------

t0 <- Sys.time()
results <- parallel::mclapply(seq_len(nrow(mutations)), compare_distributions, mc.cores = 10)
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) stop("failed mutants: ", paste(which(failed), collapse = ", "), "\n", results[[which(failed)[1]]])
cat("mutants:", length(results), " elapsed:", format(Sys.time() - t0), "\n")

stack <- function(level, name) do.call(rbind, lapply(results, function(r) r[[level]][[name]]))
metric_names <- names(results[[1]]$site)
scan <- list(
  mutations = mutations,
  lambda = get_evalue(wt),
  site = setNames(lapply(metric_names, function(name) stack("site", name)), metric_names),
  mode = setNames(lapply(metric_names, function(name) stack("mode", name)), metric_names),
  nu_range = do.call(rbind, lapply(results, `[[`, "nu_range"))
)
for (level in c("site", "mode")) for (name in metric_names) {
  stopifnot(all(is.finite(scan[[level]][[name]])))
}
saveRDS(scan, file.path(out_dir, "metrics_scan.rds"))
cat("metrics:", paste(metric_names, collapse = ", "), "\n")
