# site-by-site response matrices ------------------------------------------

#' Calculate mode-by-site motion-response matrices
#'
#' @param mutants A tibble of single-point mutants generated using `generate_mutants`
#'
#' @return a response matrix of the form \eqn{R_{nj}} (response mode is \code{n}, mutated site is \code{j})
#'
#' @name mrs_motion_by_mode
#'
NULL

# Ensemble mode-site response matrices----------------------------------------------


#' @rdname mrs_motion_by_mode
#'
#' @details  `mrs_motion_dmsfnj()` calculates the change of msf along mode n averaged over mutations at site j.
#'
#' @export
#'
mrs_motion_dmsfnj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           dmsfnjm = map2(wt, mut, delta_motion_dmsfn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, dmsfnjm)) %>%
    select(n, j, mutation, dmsfnjm) %>%
    group_by(n, j) %>%
    summarise(dmsfnj = mean(dmsfnjm)) %>%
    ungroup()

  result
}

#' @rdname mrs_motion_by_mode
#'
#' @details  `mrs_motion_dhnj()` calculates change of entropy contribution of mode n averaged over mutations at site j.
#'
#' @export
#'
mrs_motion_dhnj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           dhnjm = map2(wt, mut, delta_motion_dhn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, dhnjm)) %>%
    select(n, j, mutation, dhnjm) %>%
    group_by(n, j) %>%
    summarise(dhnj = mean(dhnjm)) %>%
    ungroup()

  result
}



#' @rdname mrs_motion_by_mode
#'
#' @details  `mrs_motion_nhnj()` calculates conservation score nh for mode n averaged over mutations at site j.
#'
#' @export
#'
mrs_motion_nhnj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           nhnjm = map2(wt, mut, delta_motion_nhn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, nhnjm)) %>%
    select(n, j, mutation, nhnjm) %>%
    group_by(n, j) %>%
    summarise(nhnj = mean(nhnjm)) %>%
    ungroup()

  result
}

#' @rdname mrs_motion_by_mode
#'
#' @details  `mrs_motion_rwsipnj()` calculates rwsip for mode n averaged over mutations at j
#'
#' @export
#'
mrs_motion_rwsipnj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           rwsipnjm = map2(wt, mut, delta_motion_rwsipn)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, rwsipnjm)) %>%
    select(n, j, mutation, rwsipnjm) %>%
    group_by(n, j) %>%
    summarise(rwsipnj = mean(rwsipnjm)) %>%
    ungroup()

  result
}






