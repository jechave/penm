#' Initialize prot object
init_prot <- function(prot,  model, d_max,  frustrated) {

  result <- prot
  result <- add_site_indexes(result)
  result <- add_enm(result, model, d_max, frustrated)
  result$energy <- enm_energy(result)

  result
}


#' add site indexes
#'
add_site_indexes <- function(prot) {
  result <- prot
  result$nsites <- length(result$pdb_site) # add nsites
  result$site <- seq(result$nsites) # add site
  result
}


