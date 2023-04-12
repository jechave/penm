#' Calculate ddgact profile using smrs method
#'
#' An influence profile is a vector with elements (Vj), where Vj is a total protein response to mutations of site j, average over mutations of j.
#' It uses a simulation method (calculates responses for various instances of forces, then calculates averages)
#' `ddgact` is the mutational change in activation free energy calculated using the lfenm model
#' It is the energy necessary to deform the mutant's active site to "push it back" to its wild-type conformation
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
#' @param pdb_site_active is a list of sites (using their resno in the pdb file) that define an "active site"
#' @param nmut is the number of mutations per site to simulate
#' @param mut_model is the mutational model, which may be "lfenm" or "sclfenm"
#' @param mut_dl_sigma is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).
#' @param mut_sd_min is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated
#' @param seed is an integer seed to set_seed before generating random mutations
#'
#' @return a tibble containing the ddg profile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pdb <- bio3d::read.pdb("2acy")
#' wt <- set_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' pdb_site_active <- c(10,14,21)
#' responses <- smrs_ddgact(wt, pdb_site_active, nmut = 10, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1024)
#' }
#'
#' @family mutscan functions
#'
smrs_ddgact <- function(wt, pdb_site_active, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  stopifnot(mut_model == "lfenm") # no ddgact_tds implementation yet

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut, mut_model, mut_dl_sigma, mut_sd_min)

  # generate perturbations

  perturbations <- generate_perturbations(wt, nmut, mut_dl_sigma, mut_sd_min, seed)

  dlmat <- perturbations$dlmat # the matrix of dlij perturbations (dlij is the change of spring length of contact i-j)
  fmat <- perturbations$fmat  # matrix of forces (fij = kij * dlij, directed along the contact)

  # ddgact
  amat <- amat_asite(wt, pdb_site_active) # matrix sqrt(t(C) Kaa C)), to pass to calculate_s2ij
  ddgact_dv  <- 1/2 * colSums(calculate_s2ij(fmat, amat))
  ddgact_tds <- 0 # assume "lfenm" mut_model
  ddgact <- ddgact_dv + ddgact_tds

  result <- tibble(
    site = get_site(wt),
    pdb_site = get_pdb_site(wt),
    is_active = pdb_site %in% pdb_site_active,
    ddgact
  )

  result
}


#' Caculate effective K matrix of active site
#'
#' @noRd
#'
amat_asite <- function(prot, pdb_site_active) {
  asite <- active_site_indexes(prot, pdb_site_active)
  cmat <- get_cmat(prot)
  cmat_asite <- cmat[asite$ind_active, asite$ind_active]
  kmat_asite <- solve(cmat_asite)

  kaa_full <- matrix(0, nrow(cmat), ncol(cmat))
  kaa_full[asite$ind_active, asite$ind_active] <-  kmat_asite

  amat2 <- t(cmat) %*% (kaa_full %*% cmat)
  amat <- matrix_sqrt(amat2)
  amat
}

#' Calculate \eqn{\mathbf{C}^{1/2}}
#'
#' Calculates the (matrix) square root of a matrix \code{m}
#'
#' @param m is a real square matrix
#' @returns the sqrt of the matrix
#'
#' @noRd
#'
#' @family mutscan helpers
matrix_sqrt <- function(m) {
  eig <- eigen(m)

  evalue <- round(eig$values, 14) # to make zero small negative numbers with double precission errors
  stopifnot(min(evalue) >= 0)

  umat <- eig$vectors

  result <- umat %*% (sqrt(evalue) * t(umat))
  result
}




