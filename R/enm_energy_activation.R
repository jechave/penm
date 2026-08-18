#' Activation free energy, internal energy contribution
#'
#' Calculates the energy cost of deforming the active site of \code{prot} into the
#' conformation it has in \code{ideal}, treating the transition state as rigid.
#'
#' @param prot A protein object obtained using \code{set_enm()}
#' @param ideal A protein object whose conformation defines the ideal (transition) state
#' @param pdb_site_active A vector of active-site residues, in pdb numbering (\code{resno}).
#'   If \code{NA} (the default), the function returns \code{NA}.
#'
#' @return A scalar, the internal-energy contribution to the activation free energy,
#'   or \code{NA} if \code{pdb_site_active} is \code{NA}.
#'
#' @export
#' @family enm_energy
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


#' Activation free energy, entropic contribution
#'
#' Calculates the entropic contribution to the activation free energy, as the
#' difference between the entropic free energy of the active site of \code{prot} and
#' that of the transition state, which is assumed rigid (contributing zero).
#'
#' @param prot A protein object obtained using \code{set_enm()}
#' @param ideal A protein object whose conformation defines the ideal (transition) state
#' @param pdb_site_active A vector of active-site residues, in pdb numbering (\code{resno}).
#'   If \code{NA} (the default), the function returns \code{NA}.
#' @param beta Inverse temperature, \code{1/kT}
#'
#' @return A scalar, the entropic contribution to the activation free energy,
#'   or \code{NA} if \code{pdb_site_active} is \code{NA}.
#'
#' @export
#' @family enm_energy
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


#' Caculate effective K matrix of active site
#'
#' @noRd
#'
kmat_asite <- function(prot, pdb_site_active) {
  asite <- active_site_indexes(prot, pdb_site_active)
  cmat <- get_cmat(prot)
  cmat_asite <- cmat[asite$ind_active, asite$ind_active]
  kmat_asite <- solve(cmat_asite)
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



