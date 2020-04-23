#' Calculate kij for the ANM
#'
kij_anm  <- function(dij, d_max = 10, k = 1,  ...) {
  ifelse(dij <= d_max, k, 0)
}

#' Calculate kij for the GNM
#'
kij_gnm <- function(dij, d_max = 10, k = 1, ...) {
  ifelse(dij <= d_max, k, 0)
}

#' Calculate kij for model by Hinsen
#'
kij_hnm <- function(dij, ...){
  ab <- 860
  b <-  2390
  al <- 1280000
  c <- 4
  if (dij <= c) {
    kij <- ab * dij - b
  } else {
    kij <- al / dij^6
  }
  kij
}

#' Calculate kij for exponential model Hinsen
#'
kij_hnm0 <- function(dij, c = 7.5, a = 1, ...) {
  a * exp( -(dij / c) ^ 2)
}

#' Calculate kij for model by Ming and Wall (2005)
#'
kij_ming_wall <- function(dij, sdij,
                          d_max = 10.5, k = 4.5, a = 42, ...) {
  kij <-  ifelse(dij <=  d_max, k, 0)
  kij[sdij == 1] <- a * kij[sdij == 1]
  kij
}

#' Calculate kij for parameter-free anm (by Yang et al.)
#'
kij_pfgnm <- function(dij, ...) {
  1/dij^2
}

#' Calculate kij for model by Reach et al.
#'
kij_reach <- function(dij, sdist = 5, same.chain = T, ...) {
  k12 <- 712
  k13 <- 6.92
  k14 <- 32.0
  ain <- 2560
  bin <- 0.8
  aex <- 1630
  bex <- 0.772
  if (sdist == 1) {
    kij <- k12
  } else if (sdist == 2) {
    kij <- k13
  } else if (sdist == 3) {
    kij <- k14
  } else if (same.chain) {
    kij <- ain * exp(-bin * dij)
  } else {
    kij <- aex * exp(-bex * dij)
  }
  kij
}
