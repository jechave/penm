#' Calculate ddgact profile using amrs method
#'
#' `ddgact` is the mutational change in activation free energy calculated using the lfenm model
#' It is the energy necessary to deform the mutant's active site to "push it back" to its wild-type conformation.
#' \code{amrs_ddgact()} calculates a profile of site-specific \code{ddgact} values, where the jth value corresponds to averaging over mutations at site j.
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
#' @param pdb_site_active is a list of sites (using their resno in the pdb file) that define an "active site"
#' @param mut_model is the mutational model, which may be "lfenm" or "sclfenm"
#' @param mut_dl_sigma is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).
#' @param mut_sd_min is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated
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
#' responses <- amrs_ddgact(wt, pdb_site_active,  mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1)
#' }
#'
#' @family mutscan functions
#'
amrs_ddgact <- function(wt, pdb_site_active, mut_model, mut_dl_sigma, mut_sd_min) {

  stopifnot(mut_model == "lfenm") # no ddgact_tds implementation yet

  # generate perturbations

  # ddgact
  amat <- amat_asite(wt, pdb_site_active) # matrix sqrt(t(C) Kaa C)), to pass to calculate_s2ij
  ddgact_dv <- 1/2 * colSums(calculate_Rij_amrs(wt, amat, mut_dl_sigma, mut_sd_min))
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

