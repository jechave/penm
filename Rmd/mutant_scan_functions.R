


# get mutant table
get_mutants_table <- function(wt, nmut_per_site) {
  mutation <- seq(from = 0, to = nmut_per_site)
  j <- get_site(wt)
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
           dr2nj = map2(wt, mut, dr2_mode),
           de2nj = map2(wt, mut, de2_mode),
           df2nj = map2(wt, mut, df2_mode)) %>%
    select(-wt, -mut) %>%
    unnest(c(mode, dr2nj, de2nj, df2nj))
}
