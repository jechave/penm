#' Calculate kij for the ANM
#'
#' @noRd
#'
kij_anm  <- function(dij, sdij, d_max = 10, k = 1,  ...) {
  kij <- ifelse(dij <= d_max | abs(sdij) == 1, k, 0)
  kij
}

#' Calculate kij for the GNM
#'
#' @noRd
#'
kij_gnm <- kij_anm

#' Calculate kij for model by Hinsen
#'
#' @noRd
#'
kij_hnm <- function(dij, ...) {
  ab <- 860
  b <-  2390
  al <- 1280000
  c <- 4
  kij <- ifelse(dij <= c,  kij <- ab * dij - b,  kij <- al / dij^6)
  kij
}

#' Calculate kij for exponential model Hinsen
#'
#' @noRd
#'
kij_hnm0 <- function(dij, c = 7.5, a = 1, ...) {
  a * exp(-(dij / c)^2)
}

#' Calculate kij for model by Ming and Wall (2005)
#'
#' @noRd
#'
kij_ming_wall <- function(dij, sdij,
                          d_max = 10.5, k = 4.5, a = 42, ...) {
  kij <-  ifelse(dij <=  d_max, k, 0)
  kij[abs(sdij) == 1] <- a * k  #warning: regardless of distance, i,i+1 contacts are "forced"
  kij
}

#' Calculate kij for parameter-free anm (by Yang et al.)
#'
#' @noRd
#'
kij_pfanm <- function(dij, ...) {
  1 / dij^2
}



#' Calculate kij for the pfanm
#'
#' @noRd
#'
kij_pfgnm <- kij_pfanm



#' Calculate kij for model by Reach et al.
#'
#' @noRd
#'
kij_reach <- function(dij, sdist = 5, same_chain = T, ...) {
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
  } else if (same_chain) {
    kij <- ain * exp(-bin * dij)
  } else {
    kij <- aex * exp(-bex * dij)
  }
  kij
}
