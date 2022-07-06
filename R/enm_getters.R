# query prot object -------------------------------------------------------

#' Get ENM parameters
#'
#' @param prot is a prot object
#' @return list of ENM parameters
#'
#' @export
#'
get_enm_param <- function(prot) prot$param

#' Get ENM node-type
#'
#' @param prot is a prot object
#' @return node-type of ENM model (either "ca" or "sc" for alpha carbon and side-chain, respectively).
#'
#' @export
#'
get_enm_node <- function(prot)  prot$param$node


#' Get ENM model type
#'
#' @param prot is a prot object
#' @return ENM model type ("anm", "hnm0", "hnm", "ming_wall", "pfanm")
#'
#' @export
#'
get_enm_model <- function(prot)  prot$param$model

#' Get ENM contact distance cut-off
#'
#' @param prot is a prot object
#' @return d_max, the cut-off used to build the ENM network
#'
#' @export
#'
get_d_max <- function(prot) prot$param$d_max

#' Get ENM parameter frustrated
#'
#' @param prot is a prot object
#' @return frustrated, the parameter that defines whether the ENM is frustrated fully relaxed
#'
#' @export
#'
get_frustrated <- function(prot)  prot$param$frustrated


#' Get protein size
#'
#' @param prot is a prot object
#' @return number of sites (ENM nodes)
#'
#' @export
#'
get_nsites <- function(prot) prot$nodes$nsites

#' Get protein sites
#'
#' @param prot is a prot object
#' @return vector of site indexes (numbered from 1 to nsites)
#'
#' @export
#'
get_site  <- function(prot) prot$nodes$site

#' Get protein sites, pdb numeration
#'
#' @param prot is a prot object
#' @return vector of site indexes (pdb numbers)
#'
#' @export
#'
get_pdb_site <- function(prot) prot$nodes$pdb_site

#' Get the X-ray b-factor
#'
#' @param prot is a prot object
#' @return a vector of X-ray B-factors for ENM nodes
#'
#' @export
#'
get_bfactor <- function(prot) prot$nodes$bfactor

#' Get xyz coordinates
#'
#' @param prot is a prot object
#' @return A 3*N vector of xyz coordinates of the N ENM nodes
#'
#' @export
#'
get_xyz <- function(prot)  prot$nodes$xyz


#' Get ENM graph
#'
#' @param prot is a prot object
#' @return graph, a tibble containing the graph representation of the ENM
#'
#' @export
#'
get_graph <- function(prot) prot$graph

#' Get the unit vectors eij
#'
#' @param prot is a prot object
#' @return a matrix of size nedges x 3, containing unit vectors eij for all edges
#'
#' @export
#'
get_eij <- function(prot) prot$eij


#' Get the ENM network matrix K
#'
#' @param prot is a prot object
#' @return kmat, the 3N x 3N network matrix (N = nsites)
#'
#' @export
#'
get_kmat <- function(prot) prot$kmat

#' Get ENM eigenvector indexes
#'
#' @param prot is a prot object
#' @return a vector of eigenvector indexes
#'
#' @export
#'
get_mode <- function(prot) prot$nma$mode

#' Get ENM eigenvalues
#'
#' @param prot is a prot object
#' @return a vector  of eigenvalues
#'
#' @export
#'
get_evalue <- function(prot) prot$nma$evalue

#' Get ENM eigenvectors
#'
#' @param prot is a prot object
#' @return a matrix of eigenvectors U of size 3 nsites x nmodes
#'
#' @export
#'
get_umat <- function(prot) prot$nma$umat

#' Get ENM covariance matrix
#'
#' @param prot is a prot object
#' @return cmat, the covariance matrix of the ENM network, of size 3 nsites x 3 nsites
#'
#' @export
#'
get_cmat <- function(prot) prot$nma$cmat
