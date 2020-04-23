# finish this: generate list of list-elements to output
complete_prot <- function(prot,pdb_active_site,param,...) {
  res <- with(prot,{
    pdb_site_active = pdb_site_active
    n_sites = length(pdb_site)
    site = seq(n_sites)
    # set active site filters
    names(site) = pdb_site_active
    site_active = site[as.character(pdb_site_active)]
    is_active_site <- site %in% site_active
    active_site_ind <- site_to_ind(site_active)
    enm <- do.call(ganm, c(list(xyz = xyz,pdb_site = pdb_site),param$enm))
    cmat_active <- enm$cmat[active_site_ind, active_site_ind]
    kmat_active <- solve(cmat_active)
    # distance matrices
    dm <- my_dm(xyz) # euclidian distance matrix, and vectors eij and matrix Eij
    sd_ij <- my_sdm(pdb_site) # sequence distance matrix
    # set up initial graph
    v0_ij = matrix(0,n_sites,n_sites)
    graph <- set_graph(enm$k_ij,l_ij = dm$d_ij, v0_ij = v0_ij)
    lst(n_sites,
        site,
        pdb_site_active,
        is_active_site,
        active_site_ind,
        cmat = enm$cmat,
        kmat = enm$kmat,
        k_ij = enm$k_ij,
        kmat_active,
        d_ij = dm$d_ij,
        e_versor_ij = dm$e_versor_ij,
        e_matrix_ij = dm$e_matrix_ij,
        sd_ij,
        graph)
  })
  c(prot,res)
}
