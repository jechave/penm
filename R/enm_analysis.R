# Calculate various protein properties


# Calculate site profiles ----------------------------------------------------

#' Calculate CN site-dependent profile
#'
#' Calculates the Contact Number (CN) of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with cn values for each site
#'
#' @export
#'
#' @family calculate site profiles
#' @noRd
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
#' @noRd
#'
get_wcn <- function(prot) wcn_xyz(get_xyz(prot))

#' Calculate MSF site-dependent profile
#'
#' Calculates the mean-square-fluctuation of each site
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with msf values for each site
#'
#' @export
#' @noRd
#'
#' @family calculate site profiles
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
#' @noRd
#'
#' @family calculate site profiles
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

#' Site-dependent ENM minimum stress energy profile
#'
#' Calculates the sum for each site of the stress energy of each of it's springs at equilibrium
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of site-dependent stress-energy values
#'
#' @export
#' @noRd
#'
#'
#' @family calculate site profiles
#'
get_stress <- function(prot) {
  g1 <- get_graph(prot)
  g2 <- g1 %>%
    select(edge, j, i, v0ij, sdij, lij, kij, dij)
  names(g2) <- names(g1)
  g <- rbind(g1, g2)

  g <- g %>%
    mutate(stress = .5 * kij * (dij - lij)^2) %>%
    group_by(i) %>%
    summarise(stress = sum(stress))  %>%
    select(stress)

  as.vector(g$stress)

}



# Get mode profiles -------------------------------------------------------

#' Calculate MSF mode-dependent profile
#'
#' Calculates the mean-square-fluctuation in the direction of each normal mode
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a vector of size nsites with msf values for each mode
#'
#' @export
#' @noRd
#'
#' @family calculate mode profiles
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
#' @noRd
#'
#' @family calculate site-by-site matrices
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
#' @noRd
#'
#' @family calculate site-by-site matrices
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
#' @noRd
#'
#' @family calculate site-by-site matrices
#'
get_reduced_kmat <- function(prot) {
  get_kmat(prot) %>%
    reduce_matrix()
}





# site by mode matrices ---------------------------------------------------


#' Calculate MSF site-dependent profile for each mode
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nmodes with the msf of each site contributed by each mode
#'
#' @export
#' @noRd
#'
#' @family calculate site-by-mode matrices
#'
get_msf_site_mode <- function(prot) {
  umat2 <- get_umat2(prot)
  msf_site_mode <- t(t(umat2) / get_evalue(prot))
  msf_site_mode
}



#' Calculate Reduced \code{umat^2}
#'
#'Calculates a matrix of size nsites x nmodes. Element umat2(i,n) is the contribution of site i to mode n (amplitude squared, added over x,y,z)
#'
#' @param prot is a protein object obtained using set_enm()
#' @returns a matrix of size nsites x nmodes with contribution of each site to each mode.
#'
#' @export
#' @noRd
#'
#' @family calculate site-by-mode matrices
#'
get_umat2 <- function(prot) {
  umat2 <- get_umat(prot)^2
  dim(umat2) <- c(3, nrow(umat2) / 3, ncol(umat2))
  umat2 <- apply(umat2, c(2, 3), sum)
  umat2
}



