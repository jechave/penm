#' Calculate perturbation-response-scanning responses
#'
prs <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, beta) {
  mutants <- get_mutants_table(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)
  dfj <- calculate_dfj(mutants, beta)
  dfij <- calculate_dfij(mutants)
  dfnj <- calculate_dfnj(mutants)


  # add site info to dfij
  dati <- tibble(i = get_site(wt), msfi = get_msf_site(wt), mlmsi = get_mlms(wt))
  datj <- tibble(j = get_site(wt), msfj = get_msf_site(wt), mlmsj = get_mlms(wt))

  dfij <- dfij %>%
    inner_join(dati) %>%
    inner_join(datj) %>%
    select(i, j, mutation, msfi, msfj, mlmsi, mlmsj, everything())


  # add mode info to dfnj
  datn <- tibble(n = get_mode(wt), msfn = get_msf_mode(wt))
  datj <- tibble(j = get_site(wt), msfj = get_msf_site(wt), mlmsj = get_mlms(wt))

  dfnj <- dfnj %>%
    inner_join(datn) %>%
    inner_join(datj) %>%
    select(n, j, mutation, msfn, msfj, mlmsj, everything())

  lst(enm_param, mut_param, dfj,  dfij,  dfnj)
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
calculate_dfj <- function(mutants, beta) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  kmat_sqrt <- get_kmat_sqrt(wt)

  dfij <- mutants %>%
    mutate(i = map(wt, get_site),
           de2ij = map2(wt, mut, calculate_de2i, kmat_sqrt = kmat_sqrt),
           dvsij = map2(wt, mut, calculate_dvsi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, de2ij, dvsij)) %>%
    mutate(dvmij = dvsij - de2ij) %>%
    select(i, j, mutation,  de2ij, dvsij, dvmij)

  result <- dfij %>%
    group_by(j, mutation) %>%
    summarise(de2j = 1/2 * sum(de2ij),
              dvsj = 1/2 * sum(dvsij),
              dvmj = 1/2 * sum(dvmij))

  result
}


#' Calculate structural mutational response, site analysis
calculate_dfij <- function(mutants) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  kmat_sqrt <- get_kmat_sqrt(wt)

  result <- mutants %>%
    mutate(i = map(wt, get_site),
           dr2ij = map2(wt, mut, calculate_dr2i),
           df2ij = map2(wt, mut, calculate_df2i),
           de2ij = map2(wt, mut, calculate_de2i, kmat_sqrt = kmat_sqrt),
           dvsij = map2(wt, mut, calculate_dvsi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dr2ij, df2ij, de2ij, dvsij)) %>%
    mutate(dvmij = dvsij - de2ij) %>%
    select(i, j, mutation, dr2ij, df2ij, de2ij, dvsij, dvmij)

  result
}

#' Calculate structural mutational response, mode analysis
calculate_dfnj <- function(mutants) {
  # structural differences, mode analysis
  result <- mutants %>%
    mutate(n = map(wt, get_mode),
           dr2nj = map2(wt, mut, calculate_dr2n),
           de2nj = map2(wt, mut, calculate_de2n),
           df2nj = map2(wt, mut, calculate_df2n))
  result <- result %>%
    select(-wt, -mut) %>%
    unnest(c(n, dr2nj, de2nj, df2nj)) %>%
    select(n, j, mutation, dr2nj, df2nj, de2nj)
  result
}
