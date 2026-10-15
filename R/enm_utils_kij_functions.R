# Units. Every kij_* returns spring constants in kcal/(mol A^2), so that every
# energy penm computes is in kcal/mol and beta_boltzmann() (in mol/kcal) applies
# to every model. Models published in kJ/(mol A^2) (hnm, reach) keep their
# published numbers in the code and convert on return with kj_to_kcal(): the
# constants below are the papers' values, not mistyped ones. Models with no
# physical scale (anm, gnm, anm_smooth, pfanm, pfgnm, hnm0) are read as
# kcal/(mol A^2), with a scale that is a convention, not a calibration.

#' Convert a spring constant from kJ/(mol A^2) to kcal/(mol A^2)
#'
#' @param k_kj spring constants in kJ/(mol A^2)
#'
#' @returns `k_kj` in kcal/(mol A^2), using 1 kcal = 4.184 kJ
#'
#' @noRd
#'
kj_to_kcal <- function(k_kj) {
  kj_per_kcal <- 4.184
  k_kj / kj_per_kcal
}


#' Calculate kij for the ANM
#'
#' `k` is in kcal/(mol A^2). Its default of 1, also ProDy's and bio3d's
#' default, is a convention with no calibration behind it.
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
#' Hinsen et al., Chem. Phys. 261, 25 (2000), fitted to an all-atom normal-mode
#' calculation for crambin. The constants are as published, in kJ/(mol A^2)
#' (taken from bio3d's `ff.calpha`, which states those units), and the result
#' is converted to kcal/(mol A^2).
#'
#' Fitted to C-alpha networks, and negative below 2390/860 = 2.78 A, closer
#' than any two alpha carbons come; `set_enm_param()` therefore allows it only
#' for C-alpha nodes.
#'
#' @noRd
#'
kij_hnm <- function(dij, ...) {
  # published constants, kJ/(mol A^2)
  ab <- 860
  b <-  2390
  al <- 1280000
  c <- 4
  kij_kj <- ifelse(dij <= c, ab * dij - b, al / dij^6)
  kj_to_kcal(kij_kj)
}

#' Calculate kij for exponential model Hinsen
#'
#' `a` is in kcal/(mol A^2), with an arbitrary scale.
#'
#' @noRd
#'
kij_hnm0 <- function(dij, c = 7.5, a = 1, ...) {
  a * exp(-(dij / c)^2)
}

#' Calculate kij for model by Ming and Wall (2005)
#'
#' The defaults are the backbone-enhanced model of Ming & Wall, Phys. Rev. Lett.
#' 95, 198103 (2005): k = 4.26 kcal/(mol A^2), a = 42 and d_max = 10.5 A,
#' fitted to a CHARMM22 all-atom harmonic model of trypsinogen.
#'
#' @noRd
#'
kij_ming_wall <- function(dij, sdij,
                          d_max = 10.5, k = 4.26, a = 42, ...) {
  kij <-  ifelse(dij <=  d_max, k, 0)
  kij[abs(sdij) == 1] <- a * k  #warning: regardless of distance, i,i+1 contacts are "forced"
  kij
}

#' Calculate kij for parameter-free anm (by Yang et al.)
#'
#' 1 / dij^2, read as kcal/(mol A^2) with an arbitrary scale.
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
#' The constants are as published in Moritsugu et al., Biophys. J. 2009
#' (PMC2726322), Table 1, solvated myoglobin dimer, fitted to MD at 300 K, in
#' kJ/(mol A^2) (the b's in 1/A); the result is converted to kcal/(mol A^2).
#'
#' @noRd
#'
kij_reach <- function(dij, sdij, same_chain = TRUE, ...) {
  # published constants, kJ/(mol A^2)
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
  kj_to_kcal(kij)
}


#' Calculate kij for the ANM with a smoothed cutoff
#'
#' As `kij_anm()`, with the step at `d_max` replaced by
#' \eqn{\frac12 [1 - \tanh((d - d_{max}) / d_{max\_width})]}, so `kij = k/2` at the
#' cutoff, and `k` falls from its full value to 0 over about `d_max +/- 2 d_max_width`.
#' i,i+1 contacts are forced to `k` regardless of distance, as in `kij_anm`.
#' `k` is in kcal/(mol A^2), with an arbitrary scale, as in `kij_anm`.
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
#' As `kij_ming_wall()`, with the step at `d_max` replaced by
#' \eqn{\frac12 [1 - \tanh((d - d_{max}) / d_{max\_width})]}, so `kij = k/2` at the
#' cutoff, and `k` falls from its full value to 0 over about `d_max +/- 2 d_max_width`.
#' i,i+1 contacts get `a * k` regardless of distance, as in `kij_ming_wall`.
#'
#' @param d_max_width width (A) of the switch from contact to no contact around
#'   `d_max`; no default
#'
#' @noRd
#'
kij_ming_wall_smooth <- function(dij, sdij, d_max = 10.5, d_max_width, k = 4.26, a = 42, ...) {
  kij <- k * 0.5 * (1 - tanh((dij - d_max) / d_max_width))
  kij[abs(sdij) == 1] <- a * k
  kij
}
