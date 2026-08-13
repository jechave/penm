# Functions to query prot object -------------------------------------------------------


#' Get various properties of a prot object
#'
#' @param prot A protein with its associated ENM model, obtained using `set_enm`
#'
#' @name get_prot_property
#'
NULL






#' @rdname get_prot_property
#'
#' @return \code{get_enm_param}: list of ENM parameters
#'
#' @export
#'
get_enm_param <- function(prot) prot$param



#' @rdname get_prot_property
#'
#'
#' @return \code{get_enm_node}: ENM node type ("ca" or "sc")
#'
#' @noRd
#'
get_enm_node <- function(prot)  prot$param$node


#' Get ENM model type
#'
#' @param prot is a prot object
#' @return ENM model type ("anm", "hnm0", "hnm", "ming_wall", "pfanm")
#'
#'
#' @noRd
#'
get_enm_model <- function(prot)  prot$param$model


#' Get ENM contact distance cut-off
#'
#' @param prot is a prot object
#' @return d_max, the cut-off used to build the ENM network
#'
#'
#' @noRd
#'
get_d_max <- function(prot) prot$param$d_max

#' Get ENM parameter frustrated
#'
#' @param prot is a prot object
#' @return frustrated, the parameter that defines whether the ENM is frustrated fully relaxed
#'
#'
#' @noRd
#'
get_frustrated <- function(prot)  prot$param$frustrated


#' @rdname get_prot_property
#'
#'
#' @return \code{get_nsites}: number of sites (ENM nodes)
#'
#' @export
#'
get_nsites <- function(prot) prot$nodes$nsites


#' @rdname get_prot_property
#'
#'
#' @return \code{get_site}: site indexes, from `1` to `nsites`
#'
#' @export
#'
get_site  <- function(prot) prot$nodes$site


#' @rdname get_prot_property
#'
#'
#' @return \code{get_pdb_site}: site indexes, pdb numeration (resno)
#'
#' @export
#'
get_pdb_site <- function(prot) prot$nodes$pdb_site

#' @rdname get_prot_property
#'
#'
#' @return \code{get_bfactor}: B-factors of the X-ray file for the ENM nodes
#'
#' @export
#'
get_bfactor <- function(prot) prot$nodes$bfactor

#' @rdname get_prot_property
#'
#'
#' @return \code{get_xyz}: A 3*N vector of xyz coordinates of the N ENM nodes
#'
#' @export
#'

get_xyz <- function(prot)  prot$nodes$xyz


#' Get ENM graph
#'
#' @param prot is a prot object
#' @return graph, a tibble containing the graph representation of the ENM
#'
#'
#' @noRd
#'
get_graph <- function(prot) prot$graph

#' Get the unit vectors eij
#'
#' @param prot is a prot object
#' @return a matrix of size nedges x 3, containing unit vectors eij for all edges
#'
#'
#' @noRd
#'
get_eij <- function(prot) prot$eij


#' Get the ENM network matrix K
#'
#' @param prot is a prot object
#' @return kmat, the 3N x 3N network matrix (N = nsites)
#'
#' @noRd
#'
get_kmat <- function(prot) prot$kmat

#' Get ENM eigenvector indexes
#'
#' @param prot is a prot object
#' @return a vector of eigenvector indexes
#'
#'
#' @noRd
#'
get_mode <- function(prot) prot$nma$mode

#' Get ENM eigenvalues
#'
#' @param prot is a prot object
#' @return a vector  of eigenvalues
#'
#'
#' @noRd
#'
get_evalue <- function(prot) prot$nma$evalue

#' Get ENM eigenvectors
#'
#' @param prot is a prot object
#' @return a matrix of eigenvectors U of size 3 nsites x nmodes
#'
#' @noRd
#'
get_umat <- function(prot) prot$nma$umat

#' Get ENM covariance matrix
#'
#' @param prot is a prot object
#' @return cmat, the covariance matrix of the ENM network, of size 3 nsites x 3 nsites
#'
#'
#' @noRd
#'
get_cmat <- function(prot) prot$nma$cmat


#' Get number of modes
#'
#' @param prot is a prot object
#' @return number of normal modes
#'
#' @noRd
#'
get_nmodes <- function(prot) max(prot$nma$mode)
