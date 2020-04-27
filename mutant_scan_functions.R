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

get_mut <-function(wt, site, mutation = 1, dl_sigma = 1, update_enm = T, add_frust = T) {
  get_mutant_site(wt,
                  site,
                  mutation,
                  seed = 2000 + site * mutation,
                  wt, wt,
                  sd_min = param$mut$sd_min,
                  dl_sigma = dl_sigma,
                  model = param$enm$model,
                  d_max = param$enm$d_max,
                  add_frust = add_frust,
                  update_enm = update_enm)
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
    mutate(dv_min = map2_dbl(wt, mut, dv_min),
           dg_entropy = map2_dbl(wt, mut, dg_entropy),
           dv_stress = map2_dbl(wt, mut, dv_stress),
           ddv_ativation = map2_dbl(wt, mut, ddv_activation),
           dg_entropy_activation = map2_dbl(wt, mut, dg_entropy_activation)) %>%
    select(-wt, -mut)
}

#' Calculate structural mutational response, site analysis
delta_structure_site <- function(mutants) {
  # structural differences, site analysis
  mutants %>%
    mutate(i = map(wt, prot_site),
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
    mutate(mode = map(wt, prot_mode),
           dr2nj = map2(wt, mut, dr2_nm),
           de2nj = map2(wt, mut, de2_nm),
           df2nj = map2(wt, mut, df2_nm)) %>%
    select(-wt, -mut) %>%
    unnest(c(mode, dr2nj, de2nj, df2nj))
}


# Protein properties ------------------------------------------------------


# profiles
prot_site <- function(prot) prot$site
prot_pdb_site <- function(prot) prot$site
prot_bfactor <- function(prot) prot$bfactor
prot_msf <- msf_site
prot_mode <- function(prot) prot$enm$mode
prot_evalue <- function(prot) prot$enm$evalue

# scalars
prot_v_min <- function(prot) prot$energy$v_min
prot_g_entropy <- function(prot) prot$energy$g_entropy
prot_v_stress <- function(prot) prot$energy$v_stress
prot_dv_activation <- function(prot) prot$energy$dv_activation
prot_g_entropy_activation <- function(prot) prot$energy$g_entropy_activation


dv_min <- function(prot1, prot2)
  prot_v_min(prot2) - prot_v_min(prot1)

dg_entropy <- function(prot1, prot2)
  prot_g_entropy(prot2) - prot_g_entropy(prot1)

dv_stress <- function(prot1, prot2)
  prot_v_stress(prot2) - prot_v_stress(prot1)

ddv_activation <- function(prot1, prot2)
  prot_dv_activation(prot2) - prot_dv_activation(prot1)

dg_entropy_activation <- function(prot1, prot2)
  prot_g_entropy_activation(prot2) - prot_g_entropy_activation(prot1)



