# Compare protein fluctuations, site by site------------------------------------

#' Compare the fluctuations of two proteins site by site
#'
#' The fluctuations of each site are a three-dimensional Gaussian, with
#' covariance \eqn{\mathbf{C}_{ii}}{C_ii}, the site's 3 x 3 block of the
#' covariance matrix. Its principal axes describe a fluctuation ellipsoid, with
#' a size, a shape and an orientation. These functions compare the ellipsoids of
#' `wt` and `mut`, site by site. All results are in nats, and none depends on
#' the temperature or the unit of length.
#'
#' The two proteins must have the same sites.
#'
#' @param wt A protein, with normal modes
#' @param mut A second protein with the same sites, with normal modes (a genm
#'   mutant gets them from [set_enm_nma()])
#'
#' @return
#' * `delta_motion_dh_site()`: a vector of size `nsites`.
#' * `delta_motion_dlcov_site()`: a tibble with one row per site and columns
#'   `pdb_site`, `dlcov`, `size`, `shape` and `orientation`.
#'
#' @seealso [delta_motion_by_mode] for comparisons along normal modes;
#'   [get_h_site()]; [delta_structure_by_site] compares structure rather than
#'   fluctuations. [get_mutant_site()] produces the `mut` argument; [set_enm()]
#'   the `wt`.
#'
#' @references V. Arsigny, P. Fillard, X. Pennec and N. Ayache (2006).
#'   Log-Euclidean metrics for fast and simple calculus on diffusion tensors.
#'   Magnetic Resonance in Medicine 56, 411-421.
#'
#' @examples
#' # Under "genm" a mutation changes the network, and so the fluctuations.
#' # A genm mutant needs its normal modes added first.
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
#'               mut_model = "genm", d_max_graph = 14)
#' mut <- set_enm_nma(get_mutant_site(wt, site_mut = 11, mutation = 1))
#'
#' dh <- delta_motion_dh_site(wt, mut)
#' length(dh)                          # one value per site
#'
#' dlcov <- delta_motion_dlcov_site(wt, mut)
#' dlcov[which.max(dlcov$dlcov), ]
#'
#' # the three parts add up, in squares, to the whole
#' all.equal(dlcov$dlcov^2,
#'           dlcov$size^2 + dlcov$shape^2 + dlcov$orientation^2)
#'
#' # the size part is the entropy change, unsigned and rescaled
#' all.equal(dlcov$size, abs(dh) / sqrt(3))
#'
#' @name delta_motion_by_site
#'
NULL

#' @rdname delta_motion_by_site
#'
#' @details `delta_motion_dh_site()` returns the change in the entropy of each
#' site's fluctuations, `get_h_site(mut) - get_h_site(wt)`:
#' \deqn{\Delta h_i = \sum_{\alpha = 1}^{3} \ln \frac{\sigma^{mut}_{i\alpha}}{\sigma^{wt}_{i\alpha}},}{
#'       Delta h_i = sum_alpha ln(sigma^mut_i,alpha / sigma^wt_i,alpha),}
#' where the \eqn{\sigma_{i\alpha}}{sigma_i,alpha} are the standard deviations
#' along the principal axes. It is positive where `mut` fluctuates more. It
#' measures the volume of the fluctuation ellipsoid only, not its shape or
#' orientation, so the proteins need not be superposed. Summed over sites it is
#' not the change in the protein's entropy: each site's entropy ignores the
#' correlations with other sites.
#'
#' @export
#'
delta_motion_dh_site <- function(wt, mut) {
  stopifnot(get_pdb_site(wt) == get_pdb_site(mut)) # no indels
  get_h_site(mut) - get_h_site(wt)
}


#' @rdname delta_motion_by_site
#'
#' @details `delta_motion_dlcov_site()` returns, for each site, the difference
#' between the log covariances of `wt` and `mut`, split into the parts due to
#' size, shape and orientation. `dlcov` is the log-Euclidean distance (Arsigny
#' et al. 2006) between the matrices
#' \eqn{\frac{1}{2} \log \mathbf{C}_{ii}}{1/2 log C_ii}, whose eigenvalues are
#' the \eqn{\ln \sigma_{i\alpha}}{ln sigma_i,alpha}, so that a factor of two in
#' a standard deviation counts the same however large it is:
#' \deqn{\mathrm{dlcov}_i = \left\| \tfrac{1}{2} \log \mathbf{C}^{mut}_{ii} - \tfrac{1}{2} \log \mathbf{C}^{wt}_{ii} \right\|_F.}{
#'       dlcov_i = || 1/2 log C^mut_ii - 1/2 log C^wt_ii ||_F.}
#' It splits exactly into three parts,
#' \deqn{\mathrm{dlcov}^2 = \mathrm{size}^2 + \mathrm{shape}^2 + \mathrm{orientation}^2.}{
#'       dlcov^2 = size^2 + shape^2 + orientation^2.}
#' Write \eqn{\ell_\alpha = \ln \sigma_\alpha}{l_alpha = ln sigma_alpha} for each
#' protein, with the axes ordered from the longest, \eqn{\mu}{mu} for their
#' mean and \eqn{s_\alpha = \ell_\alpha - \mu}{s_alpha = l_alpha - mu} for what
#' remains, which describes the shape alone. Then
#' * `size` \eqn{= \sqrt{3} \, |\mu^{mut} - \mu^{wt}|}{= sqrt(3) |mu^mut - mu^wt|},
#'   which is \eqn{|\Delta h_i| / \sqrt{3}}{|Delta h_i| / sqrt(3)};
#' * `shape` \eqn{= \sqrt{\sum_\alpha (s^{mut}_\alpha - s^{wt}_\alpha)^2}}{= sqrt(sum_alpha (s^mut_alpha - s^wt_alpha)^2)},
#'   comparing the two ellipsoids at equal size, longest axis with longest axis;
#' * `orientation`
#'   \eqn{= \sqrt{2 \sum_{\alpha\beta} s^{wt}_\alpha s^{mut}_\beta (\delta_{\alpha\beta} - o^2_{\alpha\beta})}}{= sqrt(2 sum_alpha,beta s^wt_alpha s^mut_beta (delta_alpha,beta - o_alpha,beta^2))},
#'   where \eqn{o_{\alpha\beta}}{o_alpha,beta} is the overlap of axis
#'   \eqn{\alpha}{alpha} of `wt` with axis \eqn{\beta}{beta} of `mut`. It is
#'   zero when the axes coincide. A rotation costs in proportion to how
#'   different in length the rotated axes are: rotating two axes of nearly
#'   equal length costs little, since their directions are poorly defined, and
#'   a sphere has no orientation at all.
#'
#' The orientation depends on the frame, so the two proteins must be superposed
#' ([superpose_prot()]): the function stops if they are not, and warns if they
#' are only nearly so.
#'
#' @export
#'
delta_motion_dlcov_site <- function(wt, mut) {
  stopifnot(get_pdb_site(wt) == get_pdb_site(mut)) # no indels
  check_superposed(wt, mut)

  nsites <- get_nsites(wt)
  cmat_wt <- get_cmat(wt)
  cmat_mut <- get_cmat(mut)
  dim(cmat_wt) <- c(3, nsites, 3, nsites)   # cmat[, i, , j] is the 3 x 3 block of sites i and j
  dim(cmat_mut) <- c(3, nsites, 3, nsites)

  size <- rep(NA_real_, nsites)
  shape <- rep(NA_real_, nsites)
  orientation <- rep(NA_real_, nsites)
  for (i in seq(nsites)) {
    # principal axes, longest first
    axes_wt <- eigen(cmat_wt[, i, , i], symmetric = TRUE)
    axes_mut <- eigen(cmat_mut[, i, , i], symmetric = TRUE)
    if (any(axes_wt$values <= 0) || any(axes_mut$values <= 0)) {
      stop("the covariance of site ", i, " is not positive definite")
    }

    # ln of the standard deviation along each axis
    log_sd_wt <- 0.5 * log(axes_wt$values)
    log_sd_mut <- 0.5 * log(axes_mut$values)

    # the mean is the size; what remains, which sums to zero, is the shape
    mean_wt <- mean(log_sd_wt)
    mean_mut <- mean(log_sd_mut)
    shape_wt <- log_sd_wt - mean_wt
    shape_mut <- log_sd_mut - mean_mut

    size[i] <- sqrt(3) * abs(mean_mut - mean_wt)
    shape[i] <- sqrt(sum((shape_mut - shape_wt)^2))

    # overlap_squared[alpha, beta] is the squared overlap of axis alpha of wt with axis beta of mut
    overlap_squared <- crossprod(axes_wt$vectors, axes_mut$vectors)^2
    orientation_squared <- 2 * sum(outer(shape_wt, shape_mut) * (diag(3) - overlap_squared))
    # never negative in exact arithmetic; rounding can make a zero slightly
    # negative, and only that is set to zero. Rounding in the overlaps (~1e-16)
    # comes out of the square root as ~1e-8: a site compared with itself gets
    # an orientation of that order, not exactly zero.
    rounding <- 1e-10 * (sum(shape_wt^2) + sum(shape_mut^2))
    if (orientation_squared < -rounding) {
      stop("negative orientation term at site ", i, ": ", orientation_squared)
    }
    orientation[i] <- sqrt(max(orientation_squared, 0))
  }

  dlcov <- sqrt(size^2 + shape^2 + orientation^2)
  tibble(pdb_site = get_pdb_site(wt), dlcov = dlcov, size = size,
         shape = shape, orientation = orientation)
}
