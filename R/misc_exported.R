#' Calculate's Bolzmann's beta = 1 / RT
#'
#' @param R is Boltzmann's constant per mol
#' @param T is absolute temperature in Kelving
#'
#' @returns 1 / (RT)
#'
#' @export
#'
beta_boltzmann <- function(R = 1.986e-3, T = 298) 1/(R*T)
