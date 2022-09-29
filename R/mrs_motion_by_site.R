# site-by-site response matrices ------------------------------------------

#' Calculate site-dependent motion-response matrices
#'
#' @param mutants A tibble of single-point mutants generated using `generate_mutants`
#'
#' @return a response matrix of the form \eqn{R_{ij}} (response site is \code{i}, mutated site is \code{j})
#'
#' @name mrs_motion_by_site
#'
NULL



#' @rdname mrs_motion_by_site
#'
#' @details  `mrs_motion_dmsfij()` calculates the change of msf of site i due to mutations at j, averaged over mutations at \code{j}
#'
#' @export
#'
mrs_motion_dmsfij <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dmsfijm = map2(wt, mut, delta_motion_dmsfi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dmsfijm)) %>%
    select(i, j, mutation, dmsfijm) %>%
    group_by(i, j) %>%
    summarise(dmsfij = mean(dmsfijm)) %>%  # average over mutations
    ungroup()
  result
}


#' @rdname mrs_motion_by_site
#'
#' @details  `mrs_motion_dhij()` calculates the change in entropy of site i due to mutations at j averaged over mutations at \code{j}
#'
#' @export
#'
mrs_motion_dhij <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dhijm = map2(wt, mut, delta_motion_dhi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dhijm)) %>%
    select(i, j, mutation, dhijm) %>%
    group_by(i, j) %>%
    summarise(dhij = mean(dhijm)) %>%  # average over mutations
    ungroup()
  result
}


#' @rdname mrs_motion_by_site
#'
#' @details  `mrs_motion_rwsipij()` calculates rwsip between mutant and wt site i distributions due to mutations at j, averaged over mutations at j
#'
#' @export
#'
mrs_motion_rwsipij <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           rwsipijm = map2(wt, mut, delta_motion_rwsipi)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, rwsipijm)) %>%
    select(i, j, mutation, rwsipijm) %>%
    group_by(i, j) %>%
    summarise(rwsipij = mean(rwsipijm)) %>%  # average over mutations
    ungroup()
  result
}


#' @rdname mrs_motion_by_site
#'
#' @details  `mrs_motion_dbhatij()` calculates dbhat distance between wt and mut site i distributions, averaged over mutations at j.
#'
#' @export
#'
mrs_motion_dbhatij <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dbhatijm = map2(wt, mut, delta_motion_dbhati)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dbhatijm)) %>%
    select(i, j, mutation, dbhatijm) %>%
    group_by(i, j) %>%
    summarise(dbhatij = mean(dbhatijm)) %>%  # average over mutations
    ungroup()
  result
}



