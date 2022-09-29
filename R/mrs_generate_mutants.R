#' Generate a tibble of single-point mutants
#'
#' @param wt is the wild-type protein to mutate
#' @param nmut is the number of mutations per site to introduce
#' @param mut_model is the mutational model ("lfenm" doesn't change network K matrix,  "sclfenm" changes network K matrix)
#' @param mut_dl_sigma is the \eqn{\sigma} of a normal distribution from which dlij are obtained
#' @param mut_sd_min is the minimum sequence distance of contacts to perturbate (if sd < sd_min they're not perturbed)
#' @param seed is the seed passed to `get_mutant_site()` to generate the mutations
#'
#' @return a tibble that contains \code{nsites * nmut} mutants (mutation = 0 corresponds to wt).
#'
#' @export
#' @noRd
#'
generate_mutants <- function(wt, nmut, mut_model, mut_dl_sigma, mut_sd_min, seed) {
  mutation <- seq(from = 0, to = nmut)
  j <- get_site(wt)
  # get mutants
  mutants <-  expand_grid(wt = list(wt), j, mutation)

  mutants <- mutants %>%
    mutate(mut = pmap(list(wt, j, mutation), get_mutant_site,
                      mut_model = mut_model, mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min, seed = seed))
  mutants
}
