# Wrapper


#' Calculate all response matrices and profiles, slow simulation-based method
#'
#' This method calculates many mutants, compares wt with each mutant, then averages.
#' It is slower than smrs and much slower than amrs.
#' However, this approach can be extended to calculate "dynamical" responses (bhat, etc)
#'
#' @param wt is the wild-type protein to mutate
#' @param nmut is the number of mutations per site to introduce
#' @param mut_model is the mutational model ("lfenm" doesn't change network K matrix,  "sclfenm" changes network K matrix)
#' @param mut_dl_sigma is the \eqn{\sigma} of a normal distribution from which dlij are obtained
#' @param mut_sd_min is the minimum sequence distance of contacts to perturbate (if sd < sd_min they're not perturbed)
#' @param seed is the seed passed to `get_mutant_site()` to generate the mutations
#'
#' @return a list that contains various average response matrices and profiles
#'
#' @export
#' @noRd
#'
mrs_all <- function(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  mutants <- generate_mutants(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed)

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut, mut_model, mut_dl_sigma, mut_sd_min)

  dfij <- mutants %>%
    mrs_structure_dr2ij() %>%
    inner_join(mrs_structure_df2ij(mutants)) %>%
    inner_join(mrs_structure_de2ij(mutants)) %>%
    inner_join(mrs_structure_dvsij(mutants)) %>%
    mutate(dvmij = dvsij - de2ij)

  dfj <- dfij %>%
    group_by(j) %>%
    summarise(dr2j = sum(dr2ij),
              df2j = sum(df2ij),
              de2j = 1/2 * sum(de2ij),
              dvsj = 1/2 * sum(dvsij),
              dvmj = 1/2 * sum(dvmij)) %>%
    ungroup()

  dfi <- dfij %>%
    group_by(i) %>%
    summarise(dr2i = mean(dr2ij),
              df2i = mean(df2ij),
              de2i = mean(de2ij),
              dvsi = mean(dvsij),
              dvmi = mean(dvmij)) %>%
    ungroup()

  # structural differences, mode analysis
  dfnj <- mrs_structure_df2nj(mutants) %>%
    inner_join(mrs_structure_de2nj(mutants)) %>%
    inner_join(mrs_structure_dr2nj(mutants))

  dfn <- dfnj %>%
    group_by(n) %>%
    summarise(dr2n = mean(dr2nj),
              df2n = mean(df2nj),
              de2n = mean(de2nj)) %>%
    ungroup()

  lst(dfij, dfi, dfj, dfnj, dfn, enm_param, mut_param)
}
