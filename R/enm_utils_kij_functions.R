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
#' Pairs one, two or three apart in sequence get fixed constants; all others
#' decay exponentially with distance, with different constants within and
#' between chains. Vectorised over `dij` and `sdij`, like the other kij_*.
#'
#' @noRd
#'
kij_reach <- function(dij, sdij, same_chain = TRUE, ...) {
  k12 <- 712
  k13 <- 6.92
  k14 <- 32.0
  ain <- 2560
  bin <- 0.8
  aex <- 1630
  bex <- 0.772
  stopifnot(length(sdij) == length(dij) || length(sdij) == 1)
  if (same_chain) {
    kij <- ain * exp(-bin * dij)
  } else {
    kij <- aex * exp(-bex * dij)
  }
  sdij <- rep_len(abs(sdij), length(dij))
  kij[sdij == 1] <- k12
  kij[sdij == 2] <- k13
  kij[sdij == 3] <- k14
  kij
}


#' Calculate kij for the ANM with a smoothed cutoff
#'
#' As [kij_anm()], with the step at `d_max` replaced by
#' \eqn{\frac12 [1 - \tanh((d - d_{max}) / d_{max\_width})]}, so `kij = k/2` at the
#' cutoff, and `k` falls from its full value to 0 over about `d_max +/- 2 d_max_width`.
#' i,i+1 contacts are forced to `k` regardless of distance, as in `kij_anm`.
#'
#' @param d_max_width width (A) of the switch from contact to no contact around
#'   `d_max`; no default
#'
#' @noRd
#'
kij_anm_smooth <- function(dij, sdij, d_max = 10, d_max_width, k = 1, ...) {
  kij <- k * 0.5 * (1 - tanh((dij - d_max) / d_max_width))
  kij[abs(sdij) == 1] <- k
  kij
}


#' Calculate kij for model by Ming and Wall (2005) with a smoothed cutoff
#'
#' As [kij_ming_wall()], with the step at `d_max` replaced by
#' \eqn{\frac12 [1 - \tanh((d - d_{max}) / d_{max\_width})]}, so `kij = k/2` at the
#' cutoff, and `k` falls from its full value to 0 over about `d_max +/- 2 d_max_width`.
#' i,i+1 contacts get `a * k` regardless of distance, as in `kij_ming_wall`.
#'
#' @param d_max_width width (A) of the switch from contact to no contact around
#'   `d_max`; no default
#'
#' @noRd
#'
kij_ming_wall_smooth <- function(dij, sdij, d_max = 10.5, d_max_width, k = 4.5, a = 42, ...) {
  kij <- k * 0.5 * (1 - tanh((dij - d_max) / d_max_width))
  kij[abs(sdij) == 1] <- a * k
  kij
}
