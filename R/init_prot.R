#' Initialize prot object
init_prot <- function(prot, pdb_site_active = NA,
                      model, d_max,  frustrated) {

  prot <- add_site_indexes(prot, pdb_site_active)
  prot <- add_enm(prot, model, d_max, frustrated)

  prot
}


#' add site and active-site indexes
add_site_indexes <- function(prot, pdb_site_active) {
    pdb_site <- prot$pdb_site
    nsites <- length(pdb_site)
    site <- seq(nsites)
  if(anyNA(pdb_site_active)) {
    site_active = NA
    ind_active = NA
  } else {
    site_active <- site[pdb_site %in% pdb_site_active]
    ind_active <- site_to_ind(site_active)
  }

  prot <- c(prot,
            lst( nsites = nsites,
                 site = site,
                 pdb_site_active = pdb_site_active,
                 site_active = site_active,
                 ind_active = ind_active))
  prot
}
