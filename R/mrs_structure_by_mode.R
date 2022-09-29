# Mode by site response matrices------------------------------------------------

#' Calculate mode-by-site structure-response matrices
#'
#' @param mutants A tibble of single-point mutants generated using `generate_mutants`
#'
#' @return a response matrix of the form \eqn{R_{nj}} (response mode is \code{n}, mutated site is \code{j})
#'
#' @name mrs_structure_by_mode
#'
#'
NULL

#' @rdname mrs_structure_by_mode
#'
#' @details  `mrs_structure_df2nj()` calculates the force matrix \code{f2(n, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_df2nj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           df2njm = map2(wt, mut, delta_structure_df2n)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, df2njm)) %>%
    select(n, j, mutation, df2njm) %>%
    group_by(n, j) %>%
    summarise(df2nj = mean(df2njm)) %>%
    ungroup()

  result
}


#' @rdname mrs_structure_by_mode
#'
#' @details  `mrs_structure_de2nj()` calculates the energy matrix \code{de2(n, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_de2nj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           de2njm = map2(wt, mut, delta_structure_de2n)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, de2njm)) %>%
    select(n, j, mutation, de2njm) %>%
    group_by(n, j) %>%
    summarise(de2nj = mean(de2njm)) %>%
    ungroup()

  result
}

#' @rdname mrs_structure_by_mode
#'
#' @details  `mrs_structure_dr2nj()` calculates the structural-difference matrix \code{dr2(n, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_dr2nj <- function(mutants) {
  # structural differences, mode analysis

  result <- mutants %>%
    filter(mutation != 0) %>% # mutaiton == 0 is the wt
    mutate(n = map(wt, get_mode),
           dr2njm = map2(wt, mut, delta_structure_dr2n)) %>%
    select(-wt, -mut) %>%
    unnest(c(n, dr2njm)) %>%
    select(n, j, mutation, dr2njm) %>%
    group_by(n, j) %>%
    summarise(dr2nj = mean(dr2njm)) %>%
    ungroup()

  result
}


