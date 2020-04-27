# Protein properties ------------------------------------------------------


# profiles
get_site  <- function(prot) prot$site
get_pdb_site <- function(prot) prot$site
get_bfactor <- function(prot) prot$bfactor
get_mode <- function(prot) prot$enm$mode
get_evalue <- function(prot) prot$enm$evalue

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


# scalars
get_v_min  <- function(prot) prot$energy$v_min
get_g_entropy <- function(prot) prot$energy$g_entropy
get_v_stress <- function(prot) prot$energy$v_stress
get_dv_activation <- function(prot) prot$energy$dv_activation
get_g_entropy_activation <- function(prot) prot$energy$g_entropy_activation
