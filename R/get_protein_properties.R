# Protein properties ------------------------------------------------------


# profiles
prot_site <- function(prot) prot$site
prot_pdb_site <- function(prot) prot$site
prot_bfactor <- function(prot) prot$bfactor
prot_mode <- function(prot) prot$enm$mode
prot_evalue <- function(prot) prot$enm$evalue

#' Calculate MSF profile of prot
#'
msf_site <- function(prot) {
  stopifnot(!is.null(prot$enm$cmat))
  c <- diag(prot$enm$cmat)
  nsites <- length(c) / 3
  dim(c) <- c(3, nsites)
  msf <- colSums(c)
  msf
}


# scalars
prot_v_min <- function(prot) prot$energy$v_min
prot_g_entropy <- function(prot) prot$energy$g_entropy
prot_v_stress <- function(prot) prot$energy$v_stress
prot_dv_activation <- function(prot) prot$energy$dv_activation
prot_g_entropy_activation <- function(prot) prot$energy$g_entropy_activation
