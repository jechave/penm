get_mut_energy <- function(wt, site, mutation, mut_model,
                           seed = 2000 + site * mutation, beta = beta_boltzmann(T = 298)) {
  mut <- get_mutant_site(wt,
                         site,
                         mutation,
                         mut_model,
                         mut_dl_sigma = param$mut$mut_dl_sigma,
                         mut_sd_min = param$mut$mut_sd_min,
                         seed)
 get_mutant_site()

  tibble(v_min = enm_v_min(mut), g_entropy = enm_g_entropy(mut, beta))
}

get_mut <-function(wt, site, mutation = 1, mut_model = "sclfenm", mut_dl_sigma = 1) {
  get_mutant_site(wt,
                  site,
                  mutation,
                  mut_model,
                  mut_dl_sigma = mut_dl_sigma,
                  mut_sd_min = param$mut$mut_sd_min,
                  seed = 2000 + site * mutation)
}



# get mutant table
get_mutants_table <- function(wt, nmut_per_site) {
  mutation <- seq(from = 0, to = nmut_per_site)
  j <- wt$site
  # get mutants
  mutants <-  expand_grid(wt = list(wt), j, mutation)   %>%
    mutate(mut = pmap(list(wt, j, mutation), get_mutant_site))
  mutants
}


#' Calculate energy mutational response
delta_energy <- function(mutants) {
  # energy differences
  mutants %>%
    mutate(delta_v_min = map2_dbl(wt, mut, delta_v_min),
           delta_g_entropy = map2_dbl(wt, mut, delta_g_entropy)) %>%
    select(-wt, -mut)
}

#' Calculate structural mutational response, site analysis
delta_structure_site <- function(mutants) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  mutants %>%
    mutate(i = map(wt, get_site),
           dr2ij = map2(wt, mut, dr2_site),
           de2ij = map2(wt, mut, de2_site),
           df2ij = map2(wt, mut, df2_site)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dr2ij, de2ij, df2ij))
}

#' Calculate structural mutational response, mode analysis
delta_structure_mode <- function(mutants) {
  # structural differences, mode analysis
  mutants %>%
    mutate(mode = map(wt, get_mode),
           dr2nj = map2(wt, mut, dr2_nm),
           de2nj = map2(wt, mut, de2_nm),
           df2nj = map2(wt, mut, df2_nm)) %>%
    select(-wt, -mut) %>%
    unnest(c(mode, dr2nj, de2nj, df2nj))
}
