# simulation method -------------------------------------------------------

prs.sim <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {
  mutants <- get_mutants_table(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed)

  mutants %>%
    calculate_df2ij.prs() %>%
    inner_join(calculate_dr2ij.prs(mutants)) %>%
    inner_join(calculate_de2ij.prs(mutants))
}
