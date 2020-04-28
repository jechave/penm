#' Initialize prot object
init_prot <- function(prot, pdb_site_active = NA,
                      ideal = prot,
                      model = "ming_wall", d_max = 10.5,  frustrated = FALSE) {

  prot <- add_site_indexes(prot, pdb_site_active)
  prot <- enm_add(prot, model, d_max, frustrated)
  prot$energy <- enm_energy(prot, ideal = prot)
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
