# Calculate various protein properties


# site profiles ----------------------------------------------------

#' Calculate CN site-dependent profile
#'
#' Calculates the Contact Number (CN) of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with cn values for each site
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' cn <- get_cn(wt)
#' length(cn)                      # one value per site
#' which.max(cn)                   # most buried site by contact number
#'
get_cn <- function(prot) cn_xyz(get_xyz(prot), get_d_max(prot))

#' Calculate WCN site-dependent profile
#'
#' Calculates the Weighted Contact Number (WCN) of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with wcn values for each site
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' wcn <- get_wcn(wt)
#' length(wcn)                     # one value per site
#'
#' # WCN needs no cutoff, so unlike CN it does not depend on d_max
#' cor(wcn, get_cn(wt))
#'
get_wcn <- function(prot) wcn_xyz(get_xyz(prot))


#' Calculate distance to active site
#'
#' Calculates the distance from each site to the closest active residue
#'
#' @param prot is a protein object obtained using set_enm()
#' @param pdb_site_active is a vector of pdb resno of active residues.
#'
#' @returns a vector of size nsites with dactive values for each site
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' # pdb_site_active is in PDB numbering (resno), not sequential site index;
#' # take the values from get_pdb_site() to be sure they match
#' active <- get_pdb_site(wt)[c(10, 11, 12)]
#'
#' dact <- get_dactive(wt, active)
#' range(dact)                     # 0 at the active residues themselves
#'
get_dactive <- function(prot, pdb_site_active) {
  xyz <- get_xyz(prot)
  asite <- active_site_indexes(prot, pdb_site_active)
  site_active <- asite$site_active
  result <- dactive.xyz(xyz, site_active)
  result
}





#' Calculate MSF site-dependent profile
#'
#' Calculates the mean-square-fluctuation of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with msf values for each site
#'
#' @export
#'
#' @seealso [get_msf_mode()] for the same fluctuations resolved by mode instead.
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' msf <- get_msf_site(wt)
#' which.max(msf)                  # most mobile site
#'
#' # the site profile is the diagonal of the reduced covariance matrix
#' all.equal(msf, diag(get_reduced_cmat(wt)))
#'
get_msf_site <- function(prot) {
  diag(get_reduced_cmat(prot))
}



#' Calculate MLMS site-dependent profile
#'
#' Calculates the Mean Local Mutational Stress (MLMS) profile using graph of prot object
#'
#' @param prot is a protein object obtained using set_enm()
#' @param sdij_cut An integer cutoff of sequence distance to include in calculation
#' @returns the profile of mean-local-mutational-stress (mlms) values
#'
#' @export
#'
#' @family site profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' mlms <- get_mlms(wt)
#' length(mlms)                    # one value per site
#'
#' # sdij_cut drops near-in-sequence contacts; raising it keeps fewer springs
#' sum(get_mlms(wt, sdij_cut = 5)) < sum(mlms)
#'
get_mlms <- function(prot, sdij_cut = 2) {
  g1 <- get_graph(prot)
  g2 <- g1 %>%
    select(edge, j, i, v0ij, sdij, lij, kij, dij)
  names(g2) <- names(g1)
  g <- rbind(g1, g2)

  g <- g %>%
    filter(sdij >= sdij_cut) %>%
    group_by(i) %>%
    summarise(mlms = sum(kij))  %>%
    select(mlms)

  as.vector(g$mlms)

}



# Get mode profiles -------------------------------------------------------

#' Calculate MSF mode-dependent profile
#'
#' Calculates the mean-square-fluctuation in the direction of each normal mode
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nmodes with the msf contributed by each normal mode
#'
#' @export
#'
#' @seealso [get_msf_site()] for the same fluctuations resolved by site instead.
#'
#' @family mode profiles
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' msf_n <- get_msf_mode(wt)
#' length(msf_n) == get_nmodes(wt)
#'
#' # softest modes fluctuate most: msf is 1 / eigenvalue
#' head(sort(msf_n, decreasing = TRUE))
#'
get_msf_mode <-  function(prot) 1 / get_evalue(prot)




# get site by site matrices -----------------------------------------------


#' Calculate rho matrix
#'
#' Calculates the reduced correlation matrix (size nsites x nsites, diag(rho) = 1)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites with rho(i,j) = cmat(i,j)/sqrt(cmat(i,i) * cmat(j,j))
#'
#' @export
#'
#' @family site-by-site matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' rho <- get_rho_matrix(wt)
#' dim(rho)                        # nsites x nsites
#'
#' # it is the reduced covariance matrix normalised by its diagonal,
#' # so every site is perfectly correlated with itself
#' all.equal(unname(diag(rho)), rep(1, get_nsites(wt)))
#'
get_rho_matrix <- function(prot) {
  cmat <- get_reduced_cmat(prot)
  t(cmat / sqrt(diag(cmat))) / sqrt(diag(cmat))
}

#' Calculate reduced covariance matrix
#'
#' Calculates the reduced covariance matrix (size nsites x nsites)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites with \eqn{c_{ij} = < d\mathbf{r}_i . d\mathbf{r}_j >}
#'
#' @export
#'
#' @family site-by-site matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' cmat_r <- get_reduced_cmat(wt)
#' dim(cmat_r)                     # nsites x nsites, from the 3N x 3N cmat
#'
#' # its diagonal is the site-by-site mean-square fluctuation profile
#' all.equal(diag(cmat_r), get_msf_site(wt))
#'
get_reduced_cmat <- function(prot) {
  get_cmat(prot) %>%
    reduce_matrix()
}

#' Calculate reduced ENM K matrix
#'
#' Calculates the reduced K matrix (size nsites x nsites)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nsites \eqn{K_{ij} = Tr(\mathbf{K}_{ij})}
#'
#' @export
#'
#' @family site-by-site matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' kmat_r <- get_reduced_kmat(wt)
#' dim(kmat_r)                     # nsites x nsites, from the 3N x 3N kmat
#'
#' # off-diagonal entries are non-zero only for sites in contact
#' sum(kmat_r[upper.tri(kmat_r)] != 0)
#'
get_reduced_kmat <- function(prot) {
  get_kmat(prot) %>%
    reduce_matrix()
}





# site by mode matrices ---------------------------------------------------


#' Calculate MSF site-dependent profile for each mode
#'
#' Splits the mean-square fluctuation of each site into the contribution of each
#' normal mode. Summing over modes recovers [get_msf_site()]; summing over sites
#' recovers [get_msf_mode()].
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nmodes with the msf of each site contributed by each mode
#'
#' @export
#'
#' @family site-by-mode matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' m <- get_msf_site_mode(wt)
#' dim(m)                          # nsites x nmodes
#'
#' # the two marginals are the site and mode profiles
#' all.equal(unname(rowSums(m)), unname(get_msf_site(wt)))
#' all.equal(unname(colSums(m)), unname(get_msf_mode(wt)))
#'
get_msf_site_mode <- function(prot) {
  umat2 <- get_umat2(prot)
  msf_site_mode <- t(t(umat2) / get_evalue(prot))
  msf_site_mode
}



#' Calculate Reduced \code{umat^2}
#'
#' Calculates a matrix of size nsites x nmodes. Element umat2(i,n) is the contribution of site i to mode n (amplitude squared, added over x,y,z)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nmodes with contribution of each site to each mode.
#'
#' @export
#'
#' @family site-by-mode matrices
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' u2 <- get_umat2(wt)
#' dim(u2)                         # nsites x nmodes
#'
#' # each mode is a unit vector, so its site contributions sum to 1
#' all.equal(unname(colSums(u2)), rep(1, get_nmodes(wt)))
#'
get_umat2 <- function(prot) {
  umat2 <- get_umat(prot)^2
  dim(umat2) <- c(3, nrow(umat2) / 3, ncol(umat2))
  umat2 <- apply(umat2, c(2, 3), sum)
  umat2
}







# matrix square roots -----------------------------------------------------


#' Calculate the matrix square root of the ENM K matrix
#'
#' Calculates \eqn{\mathbf{K}^{1/2}} from the eigendecomposition of the network,
#' as \eqn{\mathbf{U} \sqrt{\lambda} \mathbf{U}^T}. A general matrix-square-root
#' routine is not usable here: \code{kmat} is singular (the six rigid-body modes
#' have zero eigenvalue), so building it from the modes is what makes it defined.
#'
#' Needed to call \code{\link{delta_structure_de2i}}, which takes
#' \code{kmat_sqrt} as an argument.
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size 3 nsites x 3 nsites, the matrix square root of \code{kmat}
#'
#' @export
#'
#' @family matrix square roots
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' ks <- get_kmat_sqrt(wt)
#' dim(ks)                         # 3 nsites x 3 nsites
#'
#' # squaring it returns the network matrix
#' all.equal(ks %*% ks, as.matrix(get_kmat(wt)), check.attributes = FALSE)
#'
get_kmat_sqrt <- function(prot) {
  evalue <- get_evalue(prot)
  umat <- get_umat(prot)
  kmat_sqrt <- umat %*% (sqrt(evalue) * t(umat))
  kmat_sqrt
}

#' Calculate the matrix square root of the ENM covariance matrix
#'
#' Calculates \eqn{\mathbf{C}^{1/2}} from the eigendecomposition of the network,
#' as \eqn{\mathbf{U} \sqrt{1/\lambda} \mathbf{U}^T}. The dual of
#' \code{\link{get_kmat_sqrt}}, which uses \eqn{\sqrt{\lambda}} instead.
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size 3 nsites x 3 nsites, the matrix square root of \code{cmat}
#'
#' @export
#'
#' @family matrix square roots
#'
#' @examples
#' wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
#'               d_max = 10.5)
#'
#' cs <- get_cmat_sqrt(wt)
#' all.equal(cs %*% cs, as.matrix(get_cmat(wt)), check.attributes = FALSE)
#'
#' # cmat is the pseudo-inverse of kmat, taken over the non-rigid modes only, so
#' # this product is not the identity but the projector onto that subspace: its
#' # trace counts the 3 nsites - 6 modes that are kept.
#' sum(diag(get_kmat_sqrt(wt) %*% cs))
#'
get_cmat_sqrt <- function(prot) {
  evalue <- get_evalue(prot)
  umat <- get_umat(prot)
  cmat_sqrt <- umat %*% (sqrt(1 / evalue) * t(umat))
  cmat_sqrt
}
