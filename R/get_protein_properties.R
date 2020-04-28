# Get scalar protein properties -------------------------------------------------------

get_v_min  <- function(prot) prot$energy$v_min

get_g_entropy <- function(prot) prot$energy$g_entropy

get_v_stress <- function(prot) prot$energy$v_stress

get_dv_activation <- function(prot) prot$energy$dv_activation

get_g_entropy_activation <- function(prot) prot$energy$g_entropy_activation



# Get site profiles ----------------------------------------------------

#' Calculate various profiles for prot
#'
get_profiles_site <-  function(prot, d_max) {

  site <- get_site(prot)
  pdb_site <- get_pdb_site(prot)
  cn <- get_cn(prot, d_max)
  wcn <- get_wcn(prot)
  bfactor <- get_bfactor(prot)
  msf <- get_msf_site(prot)

  tibble(site, pdb_site, cn, wcn, bfactor, msf)
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



# Get mode profiles -------------------------------------------------------

get_mode <- function(prot) prot$enm$mode

get_evalue <- function(prot) prot$enm$evalue

get_msf_mode <-  function(prot) 1 / get_evalue(prot)



