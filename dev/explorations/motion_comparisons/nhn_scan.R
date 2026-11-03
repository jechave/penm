# Run from the package root: Rscript dev/explorations/motion_comparisons/nhn_scan.R
out_dir <- "dev/explorations/motion_comparisons"
# Complete single-point genm scan of 2acy (vignette wild type).
# Per mutant: exact nhn, and the mutant Hessian in the wt mode basis, dK = U' (K_mut - K_wt) U,
# from which a two-level perturbative estimate of nhn is computed.
suppressMessages(devtools::load_all(quiet = TRUE))
wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5, mut_model = "genm", d_max_graph = 14,
              ensemble = 1, alphabet_size = 10, mut_dl_sigma = 0.3, mut_sd_min = 2)
u_wt <- get_umat(wt)
lambda <- get_evalue(wt)
k_wt <- as.matrix(get_kmat(wt))
nmodes <- length(lambda)
alphabet <- 0:(get_enm_param(wt)$alphabet_size - 1)
mutations <- expand.grid(letter = alphabet, site = get_site(wt))
mutations <- mutations[mutations$letter != wt$nodes$sequence[mutations$site], ]

exp_entropy <- function(w) {
  p <- w / rowSums(w)
  plogp <- p * log(p)
  plogp[p == 0] <- 0
  exp(-rowSums(plogp))
}

# weights w_nm = tan^2(theta_nm), the two-level mixing of n with m; w_nn = 1
two_level_nhn <- function(dk, level) {
  gap <- outer(level, level, "-")
  theta <- 0.5 * atan2(2 * abs(dk), abs(gap))
  w <- tan(theta)^2
  diag(w) <- 1
  exp_entropy(w)
}

one_mutant <- function(row) {
  mut <- set_enm_nma(get_mutant_site(wt, mutations$site[row], mutations$letter[row]))
  dk <- crossprod(u_wt, (as.matrix(get_kmat(mut)) - k_wt) %*% u_wt)
  list(
    nhn = delta_motion_nhn(wt, mut),
    nhn_pt = two_level_nhn(dk, lambda + diag(dk)),
    dk2 = dk^2,
    dk_diag = diag(dk),
    mu = get_evalue(mut)
  )
}

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(mutations)), one_mutant, mc.cores = 10)
stopifnot(!any(vapply(res, inherits, logical(1), "try-error")))
cat("mutants:", length(res), " elapsed:", format(Sys.time() - t0), "\n")

out <- list(
  lambda = lambda,
  mutations = mutations,
  nhn = do.call(rbind, lapply(res, `[[`, "nhn")),
  nhn_pt = do.call(rbind, lapply(res, `[[`, "nhn_pt")),
  dk_diag = do.call(rbind, lapply(res, `[[`, "dk_diag")),
  mu = do.call(rbind, lapply(res, `[[`, "mu")),
  dk2_mean = Reduce(`+`, lapply(res, `[[`, "dk2")) / length(res)
)
saveRDS(out, file.path(out_dir, "nhn_scan.rds"))
cat("any non-finite nhn:", any(!is.finite(out$nhn)), " nhn range:", range(out$nhn), "\n")
