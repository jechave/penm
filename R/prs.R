#' calcualte all fast response matrices and profiles
#'
prs <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  mutants <- get_mutants_table(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed)

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)

  dfij <- mutants %>%
    calculate_dr2ij.prs() %>%
    inner_join(calculate_df2ij.prs(mutants)) %>%
    inner_join(calculate_de2ij.prs(mutants)) %>%
    inner_join(calculate_dvsij.prs(mutants)) %>%
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
  dfnj <- calculate_df2nj.prs(mutants) %>%
    inner_join(calculate_de2nj.prs(mutants)) %>%
    inner_join(calculate_dr2nj.prs(mutants))

  dfn <- dfnj %>%
    group_by(n) %>%
    summarise(dr2n = mean(dr2nj),
              df2n = mean(df2nj),
              de2n = mean(de2nj)) %>%
    ungroup()

  lst(dfij, dfi, dfj, dfnj, dfn, enm_param, mut_param)
}


#' get mutant table
#'
get_mutants_table <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {
  mutation <- seq(from = 0, to = nmut_per_site)
  j <- get_site(wt)
  # get mutants
  mutants <-  expand_grid(wt = list(wt), j, mutation)

  mutants <- mutants %>%
    mutate(mut = pmap(list(wt, j, mutation), get_mutant_site,
                      mut_model = mut_model, mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min, seed = seed))
  mutants
}

# Site-site response matrices ---------------------------------------------




calculate_dr2ij.prs <- function(mutants) {
  result <- mutants %>%
    filter(mutation > 0) %>%  # mutation == 0  is the "no-mutation" case
    mutate(i = map(wt, get_site),
           dr2ijm = map2(wt, mut, calculate_dr2i)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dr2ijm)) %>%
    select(i, j, mutation, dr2ijm) %>%
    group_by(i, j) %>%
    summarise(dr2ij = mean(dr2ijm)) %>%  # average over mutations
    ungroup()

  result
}

calculate_df2ij.prs <- function(mutants) {
  result <- mutants %>%
    filter(mutation > 0) %>%  # mutation == 0  is the "no-mutation" case
    mutate(i = map(wt, get_site),
           df2ijm = map2(wt, mut, calculate_df2i)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, df2ijm)) %>%
    select(i, j, mutation, df2ijm) %>%
    group_by(i, j) %>%
    summarise(df2ij = mean(df2ijm)) %>%  # average over mutations
    ungroup()


  result
}

calculate_de2ij.prs <- function(mutants) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  kmat_sqrt <- get_kmat_sqrt(wt)

  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           de2ijm = map2(wt, mut, calculate_de2i, kmat_sqrt = kmat_sqrt)) %>%
    select(-wt, -mut) %>%
    unnest(c(i,  de2ijm)) %>%
    select(i, j, mutation, de2ijm) %>%
    group_by(i, j) %>%
    summarise(de2ij = mean(de2ijm)) %>%  # average over mutations
    ungroup()

  result
}

calculate_dvsij.prs <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dvsijm = map2(wt, mut, calculate_dvsi.noindel)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dvsijm)) %>%
    select(i, j, mutation, dvsijm) %>%
    group_by(i, j) %>%
    summarise(dvsij = mean(dvsijm)) %>%  # average over mutations
    ungroup()
  result
}

calculate_dvmij.prs <- function(mutants) {
  # structural differences, site analysis
  de2ij <- calculate_de2ij.prs(mutants)
  dvsij <- calculate_dvsij.prs(mutants)

  result <- inner_join(de2ij, dvsij) %>%
    mutate(dvmij = dvsij - de2ij) %>%
    select(i, j, dvmij)

  result
}



# Mode-site response matrices ---------------------------------------------


#' Calculate df2nj , mode analysis
#'
calculate_df2nj.prs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           df2njm = map2(wt, mut, calculate_df2n)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, df2njm)) %>%
    select(n, j, mutation, df2njm) %>%
    group_by(n, j) %>%
    summarise(df2nj = mean(df2njm)) %>%
    ungroup()

  result
}


#' Calculate de2nj, mode analysis
#'
calculate_de2nj.prs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           de2njm = map2(wt, mut, calculate_de2n)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, de2njm)) %>%
    select(n, j, mutation, de2njm) %>%
    group_by(n, j) %>%
    summarise(de2nj = mean(de2njm)) %>%
    ungroup()

  result
}

#' Calculate dr2nj , mode analysis
#'
calculate_dr2nj.prs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           dr2njm = map2(wt, mut, calculate_dr2n)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, dr2njm)) %>%
    select(n, j, mutation, dr2njm) %>%
    group_by(n, j) %>%
    summarise(dr2nj = mean(dr2njm)) %>%
    ungroup()

  result
}



# pair comparison ---------------------------------------------------------

calculate_dvsi.noindel <- function(wt, mut) {
  gwt <- get_graph(wt)
  gmut <- get_graph(mut)

  stopifnot(all(gmut$edge == gwt$edge)) # for "lfenm": this works if the network didn't change its topology

  gmut$vsij = 1/2 * gmut$kij * (gwt$dij - gmut$lij)^2
  gwt$vsij = 1/2 * gwt$kij * (gwt$dij - gwt$lij)^2
  dvsij = gmut$vsij - gwt$vsij

  dvsij_non_zero <- !near(dvsij, 0)

  dvsij <- dvsij[dvsij_non_zero]
  i <- gwt$i[dvsij_non_zero]
  j <- gwt$j[dvsij_non_zero]
  sites_non_zero <- unique(c(i,j))

  dvsi <- rep(0, get_nsites(wt))

  for (e in seq_along(dvsij))  {
    dvsi[i[e]] = dvsi[i[e]] + dvsij[e]
    dvsi[j[e]] = dvsi[j[e]] + dvsij[e]
  }

  dvsi
}









