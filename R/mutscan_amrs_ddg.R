#' Calculate site-dependent ddg profile using amrs method
#'
#' \code{amrs_ddg()} uses an analytical method to calculate averages over mutations of elastic-network-model stability changes.
#'
#' If code{ddg[j,m]} is the mutational change in free energy calculated using the lfenm model for a mutation `m` at site `j`,
#' then this function returns a vector with values \code{ddg[j] = rowMeans(ddg[j,m])}, i.e., averages over mutations.
#'
#' For details see \doi{10.7717/peerj.11330}
#'
#' @param wt is the (wild-type) protein to mutate (an object obtained using \code{set_enm})
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
#' responses <- amrs_ddg(wt, mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1, seed = 1024)
#' }
#'
#' @family mutscan functions
#'
amrs_ddg <- function(wt, mut_model, mut_dl_sigma, mut_sd_min) {

  stopifnot(mut_model == "lfenm") # no ddg_tds implementation yet

  # ddg
  dat <-   calculate_de2ij_amrs(wt, mut_dl_sigma, mut_sd_min) %>%
    inner_join(calculate_dvsij_amrs(wt, mut_dl_sigma, mut_sd_min)) %>%
    mutate(dvmij = dvsij - de2ij)

  result <- dat %>%
    group_by(j) %>%
    summarise(
      dvsj = 1/2 * sum(dvsij),
      de2j = 1/2 * sum(de2ij)) %>%
    mutate(
      ddg_dv = dvsj - de2j,
    ) %>%
    transmute(
      site = get_site(wt),
      pdb_site = get_pdb_site(wt),
      ddg = ddg_dv
    )

  result
}

