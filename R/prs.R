#' Calculate perturbation-response-scanning responses
#'
prs <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, beta) {
  mutants <- get_mutants_table(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)
  dej <- delta_energy(mutants, beta)
  dfij <- delta_structure_site(mutants)
  dfnj <- delta_structure_mode(mutants)

  # add j-site info to dej
  dfj <- tibble(j = get_site(wt), msfj = get_msf_site(wt), mlmsj = get_mlms(wt))

  dej <- dej %>%
    inner_join(dfj)

  # add site info to dfij
  dfi <- tibble(i = get_site(wt), msfi = get_msf_site(wt), mlmsi = get_mlms(wt))
  dfj <- tibble(j = get_site(wt), msfj = get_msf_site(wt), mlmsj = get_mlms(wt))

  dfij <- dfij %>%
    inner_join(dfi) %>%
    inner_join(dfj) %>%
    select(i, j, mutation, msfi, msfj, mlmsi, mlmsj, everything())


  # add mode info to dfnj
  dn <- tibble(mode = get_mode(wt), msfn = get_msf_mode(wt))
  dj <- tibble(j = get_site(wt), msfj = get_msf_site(wt), mlmsj = get_mlms(wt))

  dfnj <- dfnj %>%
    inner_join(dn) %>%
    inner_join(dj) %>%
    select(mode, j, mutation, msfn, msfj, mlmsj, everything())

  lst(enm_param, mut_param, dej,  dfij,  dfnj)
}



# get mutant table
get_mutants_table <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min) {
  mutation <- seq(from = 0, to = nmut_per_site)
  j <- get_site(wt)
  # get mutants
  mutants <-  expand_grid(wt = list(wt), j, mutation)

  mutants <- mutants %>%
    mutate(mut = pmap(list(wt, j, mutation), get_mutant_site,
                      mut_model = mut_model, mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min))
  mutants
}


#' Calculate energy mutational response
delta_energy <- function(mutants, beta) {
  # energy differences
  mutants %>%
    mutate(dvm = map2_dbl(wt, mut, calculate_dvm),
           dvs = map2_dbl(wt, mut, calculate_dvs),
           delta_g_entropy = map2_dbl(wt, mut, delta_g_entropy, beta = beta)) %>%
    select(-wt, -mut)
}


#' Calculate structural mutational response, site analysis
delta_energy_site <- function(mutants) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  kmat_sqrt <- get_kmat_sqrt(wt)

  result <- mutants %>%
    mutate(i = map(wt, get_site),
           de2ij = map2(wt, mut, calculate_de2i, kmat_sqrt = kmat_sqrt))

  # result <- result %>%
  #   mutate(dvmij = map2(wt, mut, dvm_site))

  result <- result %>%
    mutate(dvsij = map2(wt, mut, calculate_dvsi))

  result <- result %>%
     mutate(dvmij = map2(dvsij, de2ij, ~ .x - .y))


  result %>%
    select(-wt, -mut) %>%
    unnest(c(i,  dvsij, dvmij, de2ij)) %>%
    select(i, j, mutation, de2ij, dvsij, dvmij)
}



#' Calculate structural mutational response, site analysis
delta_structure_site <- function(mutants) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  kmat_sqrt <- get_kmat_sqrt(wt)

  result <- mutants %>%
    mutate(i = map(wt, get_site),
           dr2ij = map2(wt, mut, calculate_dr2i))

  result <- result %>%
    mutate(de2ij = map2(wt, mut, calculate_de2i, kmat_sqrt = kmat_sqrt))

  result <- result %>%
    mutate(df2ij = map2(wt, mut, calculate_df2i))

  result %>%
    select(-wt, -mut) %>%
    unnest(c(i, dr2ij, de2ij, df2ij))
}

#' Calculate structural mutational response, mode analysis
delta_structure_mode <- function(mutants) {
  # structural differences, mode analysis
  result <- mutants %>%
    mutate(mode = map(wt, get_mode),
           dr2nj = map2(wt, mut, calculate_dr2n),
           de2nj = map2(wt, mut, calculate_de2n),
           df2nj = map2(wt, mut, calculate_df2n))
  result <- result %>%
    select(-wt, -mut) %>%
    unnest(c(mode, dr2nj, de2nj, df2nj))
}
