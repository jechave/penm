# query prot object -------------------------------------------------------

get_enm_node <- function(prot)  prot$param$node

get_enm_model <- function(prot)  prot$param$model

get_d_max <- function(prot) prot$param$d_max

get_frustrated <- function(prot)  prot$param$frustrated


get_nsites <- function(prot) prot$nodes$nsites

get_site  <- function(prot) prot$nodes$site

get_pdb_site <- function(prot) prot$nodes$pdb_site

get_bfactor <- function(prot) prot$nodes$bfactor

get_xyz <- function(prot)  prot$nodes$xyz


get_graph <- function(prot) prot$graph

get_eij <- function(prot) prot$eij

get_kmat <- function(prot) prot$kmat

get_mode <- function(prot) prot$nma$mode

get_evalue <- function(prot) prot$nma$evalue

get_umat <- function(prot) prot$nma$umat

get_cmat <- function(prot) prot$nma$cmat







