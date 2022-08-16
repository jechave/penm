#' Calculate all response matrices and profiles, "simulation" prs method
#'
#' @param wt is the wild-type protein to mutate
#' @param nmut_per_site is the number of mutations per site to introduce
#' @param mut_model is the mutational model ("lfenm" doesn't change network K matrix,  "sclfenm" changes network K matrix)
#' @param mut_dl_sigma is the \eqn{\sigma} of a normal distribution from which dlij are obtained
#' @param mut_sd_min is the minimum sequence distance of contacts to perturbate (if sd < sd_min they're not perturbed)
#' @param seed is the seed passed to `get_mutant_site()` to generate the mutations
#'
#' @return a list that contains various average response matrices and profiles
#'
#' @export
#'
mrs_all <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {

  mutants <- generate_mutants(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed)

  enm_param <- get_enm_param(wt)
  mut_param <- lst(nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min)

  dfij <- mutants %>%
    calculate_dr2ij_mrs() %>%
    inner_join(calculate_df2ij_mrs(mutants)) %>%
    inner_join(calculate_de2ij_mrs(mutants)) %>%
    inner_join(calculate_dvsij_mrs(mutants)) %>%
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
  dfnj <- calculate_df2nj_mrs(mutants) %>%
    inner_join(calculate_de2nj_mrs(mutants)) %>%
    inner_join(calculate_dr2nj_mrs(mutants))

  dfn <- dfnj %>%
    group_by(n) %>%
    summarise(dr2n = mean(dr2nj),
              df2n = mean(df2nj),
              de2n = mean(de2nj)) %>%
    ungroup()

  lst(dfij, dfi, dfj, dfnj, dfn, enm_param, mut_param)
}



#' Generate a tibble of single-point mutants
#'
#' @param wt is the wild-type protein to mutate
#' @param nmut_per_site is the number of mutations per site to introduce
#' @param mut_model is the mutational model ("lfenm" doesn't change network K matrix,  "sclfenm" changes network K matrix)
#' @param mut_dl_sigma is the \eqn{\sigma} of a normal distribution from which dlij are obtained
#' @param mut_sd_min is the minimum sequence distance of contacts to perturbate (if sd < sd_min they're not perturbed)
#' @param seed is the seed passed to `get_mutant_site()` to generate the mutations
#'
#' @return a tibble that contains \code{nsites * nmut_per_site} mutants (mutation = 0 corresponds to wt).
#'
#' @export
#'
generate_mutants <- function(wt, nmut_per_site, mut_model, mut_dl_sigma, mut_sd_min, seed) {
  mutation <- seq(from = 0, to = nmut_per_site)
  j <- get_site(wt)
  # get mutants
  mutants <-  expand_grid(wt = list(wt), j, mutation)

  mutants <- mutants %>%
    mutate(mut = pmap(list(wt, j, mutation), get_mutant_site,
                      mut_model = mut_model, mut_dl_sigma = mut_dl_sigma, mut_sd_min = mut_sd_min, seed = seed))
  mutants
}

# Response matrices ------------------------------------------------------------



# site-by-site response matrices ------------------------------------------

#' Calculate site-dependent structural response matrices
#'
#' @param mutants A tibble of single-point mutants generated using `generate_mutants`
#'
#' @return a response matrix of the form \eqn{R_{ij}} (response site is \code{i}, mutated site is \code{j})
#'
#' @name site_mrs_matrices
#'
NULL



#' @rdname site_mrs_matrices
#'
#' @details  `calculate_df2ij_mrs()` calculates the force matrix \code{df2(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_df2ij_mrs <- function(mutants) {
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

#' @rdname site_mrs_matrices
#'
#' @details  `calculate_de2ij_mrs()` calculates the energy-difference matrix \code{de2(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_de2ij_mrs <- function(mutants) {
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

#' @rdname site_mrs_matrices
#'
#' @details  `calculate_dr2ij_mrs()` calculates the structural difference matrix \code{dr2(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_dr2ij_mrs <- function(mutants) {
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


#' @rdname site_mrs_matrices
#'
#' @details  `calculate_dvsij_mrs()` calculates the stress-energy matrix \code{dvs(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_dvsij_mrs <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dvsijm = map2(wt, mut, calculate_dvsi_same_topology)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dvsijm)) %>%
    select(i, j, mutation, dvsijm) %>%
    group_by(i, j) %>%
    summarise(dvsij = mean(dvsijm)) %>%  # average over mutations
    ungroup()
  result
}

# Ensemble responde matrices----------------------------------------------------
#' @rdname site_mrs_matrices
#'
#' @details  `calculate_dmsfij_mrs()` calculates the change of msf of site i due to mutations at j, averaged over mutations at \code{j}
#'
#' @export
#'
calculate_dmsfij_mrs <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dmsfijm = map2(wt, mut, calculate_dmsfi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dmsfijm)) %>%
    select(i, j, mutation, dmsfijm) %>%
    group_by(i, j) %>%
    summarise(dmsfij = mean(dmsfijm)) %>%  # average over mutations
    ungroup()
  result
}


#' @rdname site_mrs_matrices
#'
#' @details  `calculate_dhij_mrs()` calculates the change in entropy of site i due to mutations at j averaged over mutations at \code{j}
#'
#' @export
#'
calculate_dhij_mrs <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dhijm = map2(wt, mut, calculate_dhi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dhijm)) %>%
    select(i, j, mutation, dhijm) %>%
    group_by(i, j) %>%
    summarise(dhij = mean(dhijm)) %>%  # average over mutations
    ungroup()
  result
}


#' @rdname site_mrs_matrices
#'
#' @details  `calculate_rwsipij_mrs()` calculates rwsip between mutant and wt site i distributions due to mutations at j, averaged over mutations at j
#'
#' @export
#'
calculate_rwsipij_mrs <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           rwsipijm = map2(wt, mut, calculate_rwsipi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, rwsipijm)) %>%
    select(i, j, mutation, rwsipijm) %>%
    group_by(i, j) %>%
    summarise(rwsipij = mean(rwsipijm)) %>%  # average over mutations
    ungroup()
  result
}


#' @rdname site_mrs_matrices
#'
#' @details  `calculate_dbhatij_mrs()` calculates dbhat distance between wt and mut site i distributions, averaged over mutations at j.
#'
#' @export
#'
calculate_dbhatij_mrs <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dbhatijm = map2(wt, mut, calculate_dbhati)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dbhatijm)) %>%
    select(i, j, mutation, dbhatijm) %>%
    group_by(i, j) %>%
    summarise(dbhatij = mean(dbhatijm)) %>%  # average over mutations
    ungroup()
  result
}



# Mode by site response matrices------------------------------------------------

#' Calculate mode-dependent structural response matrices
#'
#' @param mutants A tibble of single-point mutants generated using `generate_mutants`
#'
#' @return a response matrix of the form \eqn{R_{nj}} (response mode is \code{n}, mutated site is \code{j})
#'
#' @name mode_mrs_matrices
#'
NULL

#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_df2nj_mrs()` calculates the force matrix \code{f2(n, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_df2nj_mrs <- function(mutants) {
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


#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_de2nj_mrs()` calculates the energy matrix \code{de2(n, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_de2nj_mrs <- function(mutants) {
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

#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_dr2nj_mrs()` calculates the structural-difference matrix \code{dr2(n, j)} averaged over mutations at \code{j}
#'
#' @export
#'
calculate_dr2nj_mrs <- function(mutants) {
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


# Ensemble mode-site response matrices----------------------------------------------


#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_dmsfnj_mrs()` calculates the change of msf along mode n averaged over mutations at site j.
#'
#' @export
#'
calculate_dmsfnj_mrs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           dmsfnjm = map2(wt, mut, calculate_dmsfn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, dmsfnjm)) %>%
    select(n, j, mutation, dmsfnjm) %>%
    group_by(n, j) %>%
    summarise(dmsfnj = mean(dmsfnjm)) %>%
    ungroup()

  result
}

#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_dhnj_mrs()` calculates change of entropy contribution of mode n averaged over mutations at site j.
#'
#' @export
#'
calculate_dhnj_mrs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           dhnjm = map2(wt, mut, calculate_dhn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, dhnjm)) %>%
    select(n, j, mutation, dhnjm) %>%
    group_by(n, j) %>%
    summarise(dhnj = mean(dhnjm)) %>%
    ungroup()

  result
}



#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_nhnj_mrs()` calculates conservation score nh for mode n averaged over mutations at site j.
#'
#' @export
#'
calculate_nhnj_mrs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           nhnjm = map2(wt, mut, calculate_nhn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, nhnjm)) %>%
    select(n, j, mutation, nhnjm) %>%
    group_by(n, j) %>%
    summarise(nhnj = mean(nhnjm)) %>%
    ungroup()

  result
}

#' @rdname mode_mrs_matrices
#'
#' @details  `calculate_rwsipnj_mrs()` calculates rwsip for mode n averaged over mutations at j
#'
#' @export
#'
calculate_rwsipnj_mrs <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           rwsipnjm = map2(wt, mut, calculate_rwsipn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, rwsipnjm)) %>%
    select(n, j, mutation, rwsipnjm) %>%
    group_by(n, j) %>%
    summarise(rwsipnj = mean(rwsipnjm)) %>%
    ungroup()

  result
}






