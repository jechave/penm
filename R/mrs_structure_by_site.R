# site-by-site strucutre response matrices ------------------------------------------

#' Calculate site-by-site structure-response matrices
#'
#' @param mutants A tibble of single-point mutants generated using `generate_mutants`
#'
#' @return a response matrix of the form \eqn{R_{ij}} (response site is \code{i}, mutated site is \code{j})
#'
#' @name mrs_structure_by_site
#'
NULL



#' @rdname mrs_structure_by_site
#'
#' @details  `mrs_structure_df2ij()` calculates the force matrix \code{df2(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_df2ij <- function(mutants) {
  result <- mutants %>%
    filter(mutation > 0) %>%  # mutation == 0  is the "no-mutation" case
    mutate(i = map(wt, get_site),
           df2ijm = map2(wt, mut, delta_structure_df2i)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, df2ijm)) %>%
    select(i, j, mutation, df2ijm) %>%
    group_by(i, j) %>%
    summarise(df2ij = mean(df2ijm)) %>%  # average over mutations
    ungroup()


  result
}

#' @rdname mrs_structure_by_site
#'
#' @details  `mrs_structure_de2ij()` calculates the energy-difference matrix \code{de2(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_de2ij <- function(mutants) {
  # structural differences, site analysis
  wt <- mutants$wt[[1]]
  kmat_sqrt <- get_kmat_sqrt(wt)

  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           de2ijm = map2(wt, mut, delta_structure_de2i, kmat_sqrt = kmat_sqrt)) %>%
    select(-wt, -mut) %>%
    unnest(c(i,  de2ijm)) %>%
    select(i, j, mutation, de2ijm) %>%
    group_by(i, j) %>%
    summarise(de2ij = mean(de2ijm)) %>%  # average over mutations
    ungroup()

  result
}

#' @rdname mrs_structure_by_site
#'
#' @details  `mrs_structure_dr2ij()` calculates the structural difference matrix \code{dr2(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_dr2ij <- function(mutants) {
  result <- mutants %>%
    filter(mutation > 0) %>%  # mutation == 0  is the "no-mutation" case
    mutate(i = map(wt, get_site),
           dr2ijm = map2(wt, mut, delta_structure_dr2i)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dr2ijm)) %>%
    select(i, j, mutation, dr2ijm) %>%
    group_by(i, j) %>%
    summarise(dr2ij = mean(dr2ijm)) %>%  # average over mutations
    ungroup()

  result
}


#' @rdname mrs_structure_by_site
#'
#' @details  `mrs_structure_dvsij()` calculates the stress-energy matrix \code{dvs(i, j)} averaged over mutations at \code{j}
#'
#' @export
#'
mrs_structure_dvsij <- function(mutants) {
  # structural differences, site analysis
  result <- mutants %>%
    filter(mutation > 0) %>%
    mutate(i = map(wt, get_site),
           dvsijm = map2(wt, mut, delta_structure_dvsi_same_topology)) %>%
    select(-wt, -mut) %>%
    unnest(c(i, dvsijm)) %>%
    select(i, j, mutation, dvsijm) %>%
    group_by(i, j) %>%
    summarise(dvsij = mean(dvsijm)) %>%  # average over mutations
    ungroup()
  result
}

