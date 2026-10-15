#' Activation free energy, deformation term
#'
#' Calculates the energy cost of deforming the active site of \code{prot} into the
#' conformation it has in \code{ideal}, treating the transition state as rigid:
#' \eqn{\frac12 \Delta\mathbf{x}^T \mathbf{K}_A \Delta\mathbf{x}}{(1/2) dx' K_A dx}, with
#' \eqn{\mathbf{K}_A}{K_A} the effective stiffness of the active site. The
#' activation free energy is this term plus [dgact_tds()].
#'
#' It is not exactly the internal-energy change on activation: that is this term
#' less \eqn{3 n_A / (2\beta)}{3 n_A / (2 beta)}, with \eqn{n_A}{n_A} the number of
#' active-site nodes, which the transition state cannot fluctuate. The difference
#' cancels between two proteins, so [ddgact_dv()] is exactly the change in internal
#' energy.
#'
#' @param prot A protein object obtained using \code{set_enm()}
#' @param ideal A protein object whose conformation defines the ideal (transition) state
#' @param pdb_site_active A vector of active-site residues, in pdb numbering (\code{resno}).
#'   If \code{NA} (the default), the function returns \code{NA}.
#'
#' @return A scalar, the deformation term of the activation free energy, in
#'   kcal/mol, or \code{NA} if \code{pdb_site_active} is \code{NA}.
#'
#' @export
#' @family enm_energy
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, ensemble = 7)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1)
#'
#' # active-site residues are given in pdb numbering (resno)
#' active <- get_pdb_site(wt)[c(10, 11, 12)]
#'
#' # the cost of deforming the mutant's active site back to the wild-type shape
#' dgact_dv(mut, wt, active)
#'
#' # zero when the two conformations agree, since there is nothing to deform
#' dgact_dv(wt, wt, active)
#'
dgact_dv <- function(prot, ideal, pdb_site_active = NA) {
  if (anyNA(pdb_site_active)) {
    result <- NA
  } else {
    kmat_asite <- kmat_asite(prot, pdb_site_active)
    dxyz_asite <- dxyz_asite(prot, ideal, pdb_site_active)
    result <- .5 * my_quad_form(dxyz_asite, kmat_asite, dxyz_asite)
  }
  result
}


#' Activation free energy, fluctuation term
#'
#' Calculates the free energy of freezing the fluctuations of the active site of
#' \code{prot}, the transition state being assumed rigid:
#' \eqn{\frac{1}{2\beta} \sum_a \ln \frac{2\pi}{\beta\kappa_a}}{sum_a (1/(2 beta)) ln(2 pi / (beta kappa_a))},
#' with \eqn{\kappa_a}{kappa_a} the eigenvalues of the effective stiffness of the
#' active site. The activation free energy is this term plus [dgact_dv()].
#'
#' It is not exactly the entropic part \eqn{-T\Delta S^\ddagger}{-T dS}: that is
#' this term plus \eqn{3 n_A / (2\beta)}{3 n_A / (2 beta)}, with \eqn{n_A}{n_A} the
#' number of active-site nodes. The value also depends on the unit of length (Å).
#' Both cancel between two proteins, so [ddgact_tds()] is exactly the change in
#' \eqn{-T\Delta S^\ddagger}{-T dS}.
#'
#' @param prot A protein object obtained using \code{set_enm()}
#' @param ideal A protein object whose conformation defines the ideal (transition) state
#' @param pdb_site_active A vector of active-site residues, in pdb numbering (\code{resno}).
#'   If \code{NA} (the default), the function returns \code{NA}.
#' @param beta inverse temperature \eqn{1 / RT}{1 / RT}, in mol/kcal (see [beta_boltzmann()])
#'
#' @return A scalar, the fluctuation term of the activation free energy, in
#'   kcal/mol, or \code{NA} if \code{pdb_site_active} is \code{NA}.
#'
#' @export
#' @family enm_energy
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5, ensemble = 7)
#' mut <- get_mutant_site(wt, site_mut = 11, mutation = 1)
#'
#' active <- get_pdb_site(wt)[c(10, 11, 12)]
#'
#' dgact_tds(mut, wt, active)
#'
#' # This term is built from the eigenvalue spectrum, not from the conformations.
#' # Under "lfenm" a mutation leaves the spectrum untouched, so the two calls
#' # agree exactly and the value is the same for any ideal. Under "genm" the
#' # mutation changes the spectrum, and the two calls differ.
#' dgact_tds(wt, wt, active)
#'
dgact_tds <- function(prot, ideal, pdb_site_active = NA, beta = beta_boltzmann()) {
  # Calculate entropic contribution to dg_activation
  if (anyNA(pdb_site_active)) {
    result <- NA
  } else {
    # prot, with its active site in a "relaxed" state
    kmat_asite <- kmat_asite(prot, pdb_site_active)
    eig <- eigen(kmat_asite, symmetric = TRUE, only.values = TRUE)
    evalue <- eig$values
    gact_prot <- sum(enm_g_entropy_mode(evalue, beta))
    # prot, with its active site in the ideal conformation
    gact_ideal <- 0 # assume in the TS the active-site is rigid

    result <- gact_ideal - gact_prot
  }
  result
}



## Non-exported utility functions


#' Active site indexes
#'
#' @noRd
#'
active_site_indexes <- function(prot, pdb_site_active) {
  pdb_site_active <- pdb_site_active  # active sites in pdb numbering
  site <- get_site(prot)
  pdb_site <- get_pdb_site(prot)
  site_active <- site[pdb_site %in% pdb_site_active] # active site in penm numbering
  ind_active <- xyz_indices_site(site_active) # indices of xyz coordinates of active sites
  lst(pdb_site_active, site_active, ind_active)
}


#' Calculate effective K matrix of active site
#'
#' \eqn{[(\mathbf{K}^{+})_{AA}]^{-1}}{[(K+)_AA]^-1}, a stiffness in kcal/(mol A^2):
#' built from K+, not from the covariance, which would multiply it by beta.
#'
#' @noRd
#'
kmat_asite <- function(prot, pdb_site_active) {
  asite <- active_site_indexes(prot, pdb_site_active)
  kmat_pinv <- get_kmat_pinv(prot)
  kmat_pinv_asite <- kmat_pinv[asite$ind_active, asite$ind_active]
  kmat_asite <- solve(kmat_pinv_asite)
  kmat_asite
}




#' Activation energy, internal energy term
#'
#' @noRd
#'
dxyz_asite <- function(prot, ideal, pdb_site_active = NA) {
  if (anyNA(pdb_site_active)) {
    result <- NA
  } else {
    asite <- active_site_indexes(prot, pdb_site_active)

    dxyz <- get_xyz(prot) - get_xyz(ideal)
    site_active <- asite$site_active
    dxyz <- my_as_xyz(dxyz)
    result <- as.vector(dxyz[, site_active])
  }
  result
}



