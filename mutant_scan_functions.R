get_mut_energy <- function(prot, site, mutation, update_enm, add_frust, seed = 2000 + site * mutation) {
  mut <- get_mutant_site(wt,
                         site,
                         mutation,
                         seed,
                         wt, wt,
                         sd_min = param$mut$sd_min,
                         dl_sigma = param$mut$dl_sigma,
                         model = param$enm$model,
                         d_max = param$enm$d_max,
                         add_frust = add_frust,
                         update_enm = update_enm)
  as_tibble(mut$energy)
}
