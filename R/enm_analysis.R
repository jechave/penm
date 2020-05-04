# Get scalar protein properties -------------------------------------------------------






# Get site profiles ----------------------------------------------------

#' Calculate various profiles for prot
#'
site_profiles <-  function(prot, d_max) {
  site <- get_site(prot)
  pdb_site <- get_pdb_site(prot)
  cn <- get_cn(prot, d_max)
  wcn <- get_wcn(prot)
  bfactor <- get_bfactor(prot)
  msf <- get_msf_site(prot)
  mlms <- get_mlms(prot)
  tibble(site, pdb_site, cn, wcn, bfactor, msf, mlms)
}


get_site  <- function(prot) prot$site

get_pdb_site <- function(prot) prot$site

get_bfactor <- function(prot) prot$bfactor

get_cn <- function(prot, d_max) cn_xyz(prot$xyz, d_max)

get_wcn <- function(prot) wcn_xyz(prot$xyz)

#' Calculate MSF profile of prot
#'
get_msf_site <- function(prot) {

  stopifnot(!is.null(prot$enm$cmat))
  c <- diag(prot$enm$cmat)
  nsites <- length(c) / 3
  dim(c) <- c(3, nsites)
  msf <- colSums(c)
  msf

}

#' Calculates mean local mutational stress
#'
#' The mlms of a site is the sum over its contacts of kij
#'
get_mlms <-  function(prot) {

  result <- prot$enm$kmat %>%
    reduce_matrix() %>%
    diag()

  result

}



# Get mode profiles -------------------------------------------------------

get_mode <- function(prot) prot$enm$mode

get_evalue <- function(prot) prot$enm$evalue

get_msf_mode <-  function(prot) 1 / get_evalue(prot)



# Get matrix properties ---------------------------------------------------

#' Correlation matrix
#'
rho_matrix <- function(prot) {
  cmat <- cmat(prot)
  t(cmat / sqrt(diag(cmat))) / sqrt(diag(cmat))
}

#'  Variance-covariance matrix
#'
cmat <- function(prot) {
  prot$enm$cmat %>%
    reduce_matrix()
}

#' Network K matrix
#'
kmat <- function(prot) {
  prot$enm$kmat %>%
    reduce_matrix()
}




#' MSF of each site contributed by each mode, as tibble
#'
msf_site_mode <- function(prot) {
  prot %>%
    msf_site_mode_matrix() %>%
    matrix_to_tibble(row_name = "site", col_name = "mode", value_name = "msf")
}

#' MSF of each site contributed by each mode, as matrix
#'
msf_site_mode_matrix <- function(prot) {
  result <- prot %>%
    umat2_matrix()
  result <- t(t(result) / prot$enm$evalue)
  result
}

#' Site-reduced `umat**2` matrix, as tibble
#'
umat2 <- function(prot) {
  prot %>%
    umat2_matrix() %>%
    matrix_to_tibble(row_name = "site", col_name = "mode", value_name = "umat2")
}

#' Site-reduced `umat**2` matrix, as matrix
#'
umat2_matrix <- function(prot) {
  umat2 <- prot$enm$umat^2
  dim(umat2) <- c(3, nrow(umat2) / 3, ncol(umat2))
  umat2 <- apply(umat2, c(2, 3), sum)
  umat2
}

#' turn any matrix into a tibble
matrix_to_tibble <- function(m, row_name = "i", col_name = "j", value_name = "mij") {
  result <- m %>%
    as.data.frame() %>%
    mutate(i = seq(nrow(m))) %>%
    pivot_longer(cols = 1:ncol(m),
                 names_to = "j",
                 names_prefix = "V",
                 values_to = "mij") %>%
    mutate(j = as.integer(j))  %>%
    as_tibble()

  names(result) <-  c(row_name, col_name, value_name)

  result
}








