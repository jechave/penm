#' Calculate ddg profile using smrs method
#'
#' An influence profile is a vector with elements (Vj), where Vj is a total protein response to mutations of site j, average over mutations of j.
#' It uses a simulation method (calculates responses for various instances of forces, then calculates averages)
#' `ddg` is the mutational change in free energy calculated using the lfenm model
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
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
#' responses <- smrs_ddg(wt, nmut = 10, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1024)
#' }
#'
#' @family mutscan functions
#'
smrs_ddg <- function(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  stopifnot(mut_model == "lfenm") # no ddg_tds implementation yet

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut, mut_model, mut_dl_sigma, mut_sd_min)

  # generate perturbations

  perturbations <- generate_perturbations(wt, nmut, mut_dl_sigma, mut_sd_min, seed)

  dlmat <- perturbations$dlmat # the matrix of dlij perturbations (dlij is the change of spring length of contact i-j)
  fmat <- perturbations$fmat  # matrix of forces (fij = kij * dlij, directed along the contact)


  # ddg

  de2ij <- calculate_s2ij(fmat, get_cmat_sqrt(wt))
  dvsij <- calculate_dvsij_smrs(wt, dlmat) # ojo, esto me devuelve un tibble
  dvsij <- matrix(dvsij$dvsij, get_nsites(wt), get_nsites(wt))
  dvmij <- dvsij - de2ij

  dvmj <- 1/2 * colSums(dvmij)
  ddg_dv <- dvmj
  ddg_tds <- 0 # assume "lfenm" mut_model
  ddg <- ddg_dv + ddg_tds

  result <- tibble(
    site = get_site(wt),
    pdb_site = get_pdb_site(wt),
    ddg
  )

  result
}

