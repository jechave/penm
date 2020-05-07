# Calculate various protein properties



# Get site profiles ----------------------------------------------------


get_cn <- function(prot) cn_xyz(get_xyz(prot), get_d_max(prot))

get_wcn <- function(prot) wcn_xyz(get_xyz(prot))

#' Calculate MSF profile of prot
#'
get_msf_site <- function(prot) {
  diag(get_reduced_cmat(prot))
}

#' Calculates mean local mutational stress
#'
#' The mlms of a site is the sum over its contacts of kij
#'
get_mlms <-  function(prot) {
  diag(get_reduced_kmat(prot))
}



# Get mode profiles -------------------------------------------------------


get_msf_mode <-  function(prot) 1 / get_evalue(prot)




# get site by site matrices -----------------------------------------------



#' Correlation matrix
#'
get_rho_matrix <- function(prot) {
  cmat <- get_reduced_cmat(prot)
  t(cmat / sqrt(diag(cmat))) / sqrt(diag(cmat))
}

#'  Variance-covariance matrix
#'
get_reduced_cmat <- function(prot) {
  get_cmat(prot) %>%
    reduce_matrix()
}

#' Reduced network K matrix
#'
get_reduced_kmat <- function(prot) {
  get_kmat(prot) %>%
    reduce_matrix()
}




# site by mode matrices ---------------------------------------------------



#' MSF of each site contributed by each mode, as matrix
#'
get_msf_site_mode <- function(prot) {
  result <- prot %>%
    get_umat2()
  result <- t(t(result) / get_evalue(prot))
  result
}



#' Site-reduced `umat**2` matrix, as matrix
#'
get_umat2 <- function(prot) {
  umat2 <- get_umat(prot)^2
  dim(umat2) <- c(3, nrow(umat2) / 3, ncol(umat2))
  umat2 <- apply(umat2, c(2, 3), sum)
  umat2
}










