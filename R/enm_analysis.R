# Calculate various protein properties


# site profiles ----------------------------------------------------

#' Calculate CN site-dependent profile
#'
#' Calculates the Contact Number (CN) of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with cn values for each site
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' cn <- get_cn(wt)
#' length(cn)                      # one value per site
#' which.max(cn)                   # most buried site by contact number
#'
get_cn <- function(prot) cn_xyz(get_xyz(prot), get_d_max(prot))

#' Calculate WCN site-dependent profile
#'
#' Calculates the Weighted Contact Number (WCN) of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with wcn values for each site
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' wcn <- get_wcn(wt)
#' length(wcn)                     # one value per site
#'
#' # WCN needs no cutoff, so unlike CN it does not depend on d_max
#' cor(wcn, get_cn(wt))
#'
get_wcn <- function(prot) wcn_xyz(get_xyz(prot))


#' Calculate distance to active site
#'
#' Calculates the distance from each site to the closest active residue
#'
#' @param prot is a protein object obtained using set_enm()
#' @param pdb_site_active is a vector of pdb resno of active residues.
#'
#' @returns a vector of size nsites with dactive values for each site
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' # pdb_site_active is in PDB numbering (resno), not sequential site index;
#' # take the values from get_pdb_site() to be sure they match
#' active <- get_pdb_site(wt)[c(10, 11, 12)]
#'
#' dact <- get_dactive(wt, active)
#' range(dact)                     # 0 at the active residues themselves
#'
get_dactive <- function(prot, pdb_site_active) {
  xyz <- get_xyz(prot)
  asite <- active_site_indexes(prot, pdb_site_active)
  site_active <- asite$site_active
  result <- dactive.xyz(xyz, site_active)
  result
}





#' Calculate MSF site-dependent profile
#'
#' Calculates the mean-square-fluctuation of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a vector of size nsites with the mean-square fluctuation of each
#'   site at inverse temperature `beta`, in Å²
#'
#' @export
#'
#' @seealso [get_msf_mode()] for the same fluctuations resolved by mode instead.
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' msf <- get_msf_site(wt)
#' which.max(msf)                  # most mobile site
#'
#' # the site profile is the diagonal of the reduced covariance matrix
#' all.equal(msf, diag(get_reduced_cmat(wt)))
#'
get_msf_site <- function(prot, beta = beta_boltzmann()) {
  diag(get_reduced_cmat(prot, beta))
}



#' Calculate MLMS site-dependent profile
#'
#' Calculates the Mean Local Mutational Stress (MLMS) profile using graph of prot object
#'
#' @param prot is a protein object obtained using set_enm()
#' @param sdij_cut An integer cutoff of sequence distance to include in calculation
#' @returns the profile of mean-local-mutational-stress (mlms) values: the sum of
#'   the spring constants of each site's springs, in kcal/(mol Å²)
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' mlms <- get_mlms(wt)
#' length(mlms)                    # one value per site
#'
#' # sdij_cut drops near-in-sequence contacts; raising it keeps fewer springs
#' sum(get_mlms(wt, sdij_cut = 5)) < sum(mlms)
#'
get_mlms <- function(prot, sdij_cut = 2) {
  g1 <- get_graph(prot)
  # each edge once more, from j's side
  g2 <- g1
  g2$i <- g1$j
  g2$j <- g1$i
  g <- rbind(g1, g2)

  g <- g %>%
    filter(sdij >= sdij_cut) %>%
    group_by(i) %>%
    summarise(mlms = sum(kij))  %>%
    select(mlms)

  as.vector(g$mlms)

}



# Get mode profiles -------------------------------------------------------

#' Calculate MSF mode-dependent profile
#'
#' Calculates the mean-square-fluctuation in the direction of each normal mode
#'
#' @param prot is a protein object obtained using set_enm()
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a vector of size nmodes with the mean-square fluctuation along each
#'   normal mode at inverse temperature `beta`, \eqn{1 / (\beta \lambda_n)}{1 / (beta lambda_n)}, in Å²
#'
#' @export
#'
#' @seealso [get_msf_site()] for the same fluctuations resolved by site instead.
#'
#' @family mode profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' msf_n <- get_msf_mode(wt)
#' length(msf_n) == get_nmodes(wt)
#'
#' # softest modes fluctuate most: msf is 1 / (beta * eigenvalue)
#' head(sort(msf_n, decreasing = TRUE))
#'
get_msf_mode <- function(prot, beta = beta_boltzmann()) {
  1 / (beta * get_evalue(prot))
}



#' Entropy of the fluctuations of each site
#'
#' The fluctuations of site i are a three-dimensional Gaussian, with covariance
#' \eqn{\mathbf{C}_{ii}}{C_ii}, the site's 3 x 3 block of the covariance matrix.
#' `get_h_site()` returns its differential entropy, in nats,
#' \deqn{h_i = \frac{1}{2} \ln \det (2 \pi e \, \mathbf{C}_{ii})
#'           = \frac{3}{2} \ln (2 \pi e) + \sum_{\alpha = 1}^{3} \ln \sigma_{i\alpha},}{
#'       h_i = 1/2 ln det(2 pi e C_ii) = 3/2 ln(2 pi e) + sum_alpha ln sigma_i,alpha,}
#' where the \eqn{\sigma_{i\alpha}}{sigma_i,alpha} are the standard deviations
#' along the three principal axes of the site's fluctuations. It grows with the
#' volume of the fluctuation ellipsoid, and does not depend on its orientation.
#'
#' This is the entropy of the positions of the site, not a sequence entropy.
#' Its value depends on the unit of length (Å) and on the temperature; the
#' difference between two proteins at the same temperature depends on neither
#' (see [delta_motion_dh_site()]). The \eqn{h_i}{h_i} of all sites do not add up
#' to the entropy of the protein, because they ignore the correlations between
#' sites.
#'
#' @param prot a `prot`, with normal modes
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a vector of size nsites, the entropy of each site's fluctuations, in nats
#'
#' @export
#'
#' @seealso [get_h_mode()] for the entropy along each normal mode;
#'   [delta_motion_by_site] to compare two proteins.
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' h <- get_h_site(wt)
#' which.max(h)                    # the site with the largest fluctuation volume
#'
get_h_site <- function(prot, beta = beta_boltzmann()) {
  nsites <- get_nsites(prot)
  cmat <- get_cmat(prot, beta)
  dim(cmat) <- c(3, nsites, 3, nsites)   # cmat[, i, , j] is the 3 x 3 block of sites i and j
  h <- rep(NA_real_, nsites)
  for (i in seq(nsites)) {
    block <- cmat[, i, , i]
    if (det(block) <= 0) stop("the covariance of site ", i, " is not positive definite")
    # entropy of a three-dimensional Gaussian
    h[i] <- 0.5 * log(det(2 * pi * exp(1) * block))
  }
  h
}



#' Entropy of the fluctuations along each normal mode
#'
#' The fluctuation of a protein along a direction \eqn{\mathbf{u}}{u} is a
#' one-dimensional Gaussian, with variance
#' \eqn{\sigma^2 = \mathbf{u}^T \mathbf{C} \mathbf{u}}{sigma^2 = u' C u}.
#' `get_h_mode()` returns its differential entropy, in nats,
#' \deqn{h_n = \frac{1}{2} \ln (2 \pi e \, \sigma_n^2),}{h_n = 1/2 ln(2 pi e sigma_n^2),}
#' along each normal mode \eqn{\mathbf{u}_n}{u_n} of `reference`. By default
#' the reference is the protein itself, and \eqn{\sigma_n^2}{sigma_n^2} is the
#' variance of its own mode n, [get_msf_mode()].
#'
#' With another protein as reference, the directions are the reference's modes,
#' so several proteins can be described along the same set of directions. The
#' two must then be superposed ([superpose_prot()]): the function stops if they
#' are not, and warns if they are only nearly so.
#'
#' If the two structures differ, a little of each reference mode lies along the
#' rotations of `prot`'s structure, along which `prot`'s covariance is zero. The
#' variance is taken along the part of the mode inside `prot`'s internal space:
#' \eqn{\sigma_n^2 = \mathbf{u}_n^T \mathbf{C} \mathbf{u}_n / w_n}{sigma_n^2 = u_n' C u_n / w_n},
#' where \eqn{w_n}{w_n} is the share of the mode inside that space. For
#' mutants, \eqn{1 - w_n}{1 - w_n} is tiny; the function warns if some
#' \eqn{w_n < 0.99}{w_n < 0.99} and stops if some \eqn{w_n < 0.9}{w_n < 0.9}.
#'
#' The value depends on the unit of length (Å) and on the temperature;
#' differences between proteins at the same temperature depend on neither (see
#' [delta_motion_dh_mode()]).
#'
#' @param prot a `prot`, with normal modes
#' @param reference a `prot` with normal modes and the same sites as `prot`,
#'   whose modes are the directions; by default `prot` itself
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a vector of size nmodes, the entropy of the fluctuation of `prot`
#'   along each mode of `reference`, in nats
#'
#' @export
#'
#' @seealso [get_h_site()] for the entropy of each site's fluctuations;
#'   [delta_motion_by_mode] to compare two proteins.
#'
#' @family mode profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall", d_max = 10.5,
#'               mut_model = "genm", d_max_graph = 14)
#'
#' # along its own modes, the entropy follows from the variances
#' all.equal(get_h_mode(wt), 0.5 * log(2 * pi * exp(1) * get_msf_mode(wt)))
#'
#' # a mutant, along the wild type's modes
#' mut <- set_enm_nma(get_mutant_site(wt, site_mut = 11, mutation = 1))
#' h_mut <- get_h_mode(mut, reference = wt)
#'
get_h_mode <- function(prot, reference = prot, beta = beta_boltzmann()) {
  stopifnot(get_pdb_site(reference) == get_pdb_site(prot)) # no indels
  check_superposed(reference, prot)

  directions <- get_umat(reference)
  cmat <- get_cmat(prot, beta)
  # variance of prot along each direction: the diagonal of t(U) C U
  variance <- diag(t(directions) %*% cmat %*% directions)

  # the part of a direction outside prot's internal space adds nothing to the
  # variance; dividing by the share inside gives the variance along the part inside
  share <- rowSums(squared_overlap(reference, prot))
  variance_inside <- variance / share

  # entropy of a one-dimensional Gaussian
  0.5 * log(2 * pi * exp(1) * variance_inside)
}




# get site by site matrices -----------------------------------------------


#' Calculate rho matrix
#'
#' Calculates the reduced correlation matrix (size nsites x nsites, diag(rho) = 1)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites with rho(i,j) = cmat(i,j)/sqrt(cmat(i,i) * cmat(j,j))
#'
#' @export
#'
#' @family site-by-site matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' rho <- get_rho_matrix(wt)
#' dim(rho)                        # nsites x nsites
#'
#' # it is the reduced covariance matrix normalised by its diagonal,
#' # so every site is perfectly correlated with itself
#' all.equal(unname(diag(rho)), rep(1, get_nsites(wt)))
#'
get_rho_matrix <- function(prot) {
  cmat <- get_reduced_cmat(prot)
  t(cmat / sqrt(diag(cmat))) / sqrt(diag(cmat))
}

#' Calculate reduced covariance matrix
#'
#' Calculates the reduced covariance matrix (size nsites x nsites)
#'
#' @param prot is a protein object obtained using set_enm()
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a matrix of size nsites x nsites with \eqn{c_{ij} = < d\mathbf{r}_i . d\mathbf{r}_j >}{c_ij = < dr_i . dr_j >}
#'   at inverse temperature `beta`, in Å²
#'
#' @export
#'
#' @family site-by-site matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' cmat_r <- get_reduced_cmat(wt)
#' dim(cmat_r)                     # nsites x nsites, from the 3N x 3N cmat
#'
#' # its diagonal is the site-by-site mean-square fluctuation profile
#' all.equal(diag(cmat_r), get_msf_site(wt))
#'
get_reduced_cmat <- function(prot, beta = beta_boltzmann()) {
  get_cmat(prot, beta) %>%
    reduce_matrix()
}

#' Calculate reduced ENM K matrix
#'
#' Calculates the reduced K matrix (size nsites x nsites)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites \eqn{K_{ij} = Tr(\mathbf{K}_{ij})}{K_ij = Tr(K_ij)}
#'
#' @export
#'
#' @family site-by-site matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' kmat_r <- get_reduced_kmat(wt)
#' dim(kmat_r)                     # nsites x nsites, from the 3N x 3N kmat
#'
#' # off-diagonal entries are non-zero only for sites in contact
#' sum(kmat_r[upper.tri(kmat_r)] != 0)
#'
get_reduced_kmat <- function(prot) {
  get_kmat(prot) %>%
    reduce_matrix()
}





# site by mode matrices ---------------------------------------------------


#' Calculate MSF site-dependent profile for each mode
#'
#' Splits the mean-square fluctuation of each site into the contribution of each
#' normal mode. Summing over modes recovers [get_msf_site()]; summing over sites
#' recovers [get_msf_mode()].
#'
#' @param prot is a protein object obtained using set_enm()
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a matrix of size nsites x nmodes with the msf of each site contributed
#'   by each mode, at inverse temperature `beta`, in Å²
#'
#' @export
#'
#' @family site-by-mode matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' m <- get_msf_site_mode(wt)
#' dim(m)                          # nsites x nmodes
#'
#' # the two marginals are the site and mode profiles
#' all.equal(unname(rowSums(m)), unname(get_msf_site(wt)))
#' all.equal(unname(colSums(m)), unname(get_msf_mode(wt)))
#'
get_msf_site_mode <- function(prot, beta = beta_boltzmann()) {
  umat2 <- get_umat2(prot)
  msf_site_mode <- t(t(umat2) / (beta * get_evalue(prot)))
  msf_site_mode
}



#' Calculate Reduced \code{umat^2}
#'
#' Calculates a matrix of size nsites x nmodes. Element umat2(i,n) is the contribution of site i to mode n (amplitude squared, added over x,y,z)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nmodes with contribution of each site to each mode.
#'
#' @export
#'
#' @family site-by-mode matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' u2 <- get_umat2(wt)
#' dim(u2)                         # nsites x nmodes
#'
#' # each mode is a unit vector, so its site contributions sum to 1
#' all.equal(unname(colSums(u2)), rep(1, get_nmodes(wt)))
#'
get_umat2 <- function(prot) {
  umat2 <- get_umat(prot)^2
  dim(umat2) <- c(3, nrow(umat2) / 3, ncol(umat2))
  umat2 <- apply(umat2, c(2, 3), sum)
  umat2
}







# matrix square roots -----------------------------------------------------


#' Calculate the matrix square root of the ENM K matrix
#'
#' Calculates \eqn{\mathbf{K}^{1/2}}{K^(1/2)} from the eigendecomposition of the network,
#' as \eqn{\mathbf{U} \sqrt{\lambda} \mathbf{U}^T}{U sqrt(lambda) U^T}. A general matrix-square-root
#' routine is not usable here: \code{kmat} is singular (the six rigid-body modes
#' have zero eigenvalue), so building it from the modes is what makes it defined.
#'
#' Needed to call \code{\link{delta_structure_de2i}}, which takes
#' \code{kmat_sqrt} as an argument.
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size 3 nsites x 3 nsites, the matrix square root of \code{kmat}
#'
#' @export
#'
#' @family matrix square roots
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' ks <- get_kmat_sqrt(wt)
#' dim(ks)                         # 3 nsites x 3 nsites
#'
#' # squaring it returns the network matrix
#' all.equal(ks %*% ks, as.matrix(get_kmat(wt)), check.attributes = FALSE)
#'
get_kmat_sqrt <- function(prot) {
  evalue <- get_evalue(prot)
  umat <- get_umat(prot)
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  kmat_sqrt
}

#' Calculate the matrix square root of the ENM covariance matrix
#'
#' Calculates \eqn{\mathbf{C}^{1/2}}{C^(1/2)} from the eigendecomposition of the network,
#' as \eqn{\mathbf{U} \sqrt{1/(\beta\lambda)} \mathbf{U}^T}{U sqrt(1/(beta lambda)) U^T}. The dual of
#' \code{\link{get_kmat_sqrt}}, which uses \eqn{\sqrt{\lambda}}{sqrt(lambda)} instead.
#'
#' @param prot is a protein object obtained using set_enm()
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#' @returns a matrix of size 3 nsites x 3 nsites, the matrix square root of
#'   \code{get_cmat(prot, beta)}, in Å
#'
#' @export
#'
#' @family matrix square roots
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' cs <- get_cmat_sqrt(wt)
#' all.equal(cs %*% cs, as.matrix(get_cmat(wt)), check.attributes = FALSE)
#'
#' # with beta = 1, cmat is the pseudo-inverse of kmat, taken over the non-rigid
#' # modes only, so this product is not the identity but the projector onto that
#' # subspace: its trace counts the 3 nsites - 6 modes that are kept.
#' sum(diag(get_kmat_sqrt(wt) %*% get_cmat_sqrt(wt, beta = 1)))
#'
get_cmat_sqrt <- function(prot, beta = beta_boltzmann()) {
  evalue <- get_evalue(prot)
  umat <- get_umat(prot)
  cmat_sqrt <- umat %*% (sqrt(1 / (beta * evalue)) * t(umat))
  cmat_sqrt
}
