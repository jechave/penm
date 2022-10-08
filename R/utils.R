#' Turn an xyz vector into a 3 x N matrix
#'
#' @param r is a vector of length 3 * nsites
#' @returns a matrix of dimensions 3 x nsites with the vector's elments
#'
#' @noRd
#'
my_as_xyz <- function(r) {
  r <- as.vector(r)
  n3 <- length(r)
  if(n3 %% 3 == 0) {
    n <- as.integer(length(r)/3)
    dim(r) <- c(3,n)
    r
  } else {
    print("Error in my_as_xyz: length(r) not multiple of 3")
    print(paste("length(r):",n3))
    stop("Stopping execution")
  }
}


#' Calculate wcn for all sites of protein
#'
#' @param xyz is the vector or matrix containing xyz coordinates
#'
#' @returns a vector of wcn values for all protein sites
#'
#' @noRd
#'
wcn_xyz <- function(xyz, ...) {
  xyz <- my_as_xyz(xyz)
  nsites <- ncol(xyz)
  wcn <- rep(NA,nsites)
  for (i in seq(nsites)) {
    dijv <- xyz - xyz[,i]
    dij2 <- colSums(dijv^2)
    wcn[[i]] <- sum(1 / dij2[-i])
  }
  wcn
}

#' Calculates number of contacts.
#'
#' @param xyz  xyz coordinates of sites (e.g. alpha carbons, or center of mass, etc.)
#' @param d_max cutoff radius to define contacts
#'
#' @return vector of contact numbers
#' @noRd
#'
cn_xyz <-  function(xyz, d_max) {
  xyz <- my_as_xyz(xyz)
  nsites <- ncol(xyz)
  cn <- rep(NA,nsites)
  for (i in seq(nsites)) {
    dijv <- xyz - xyz[,i]
    dij <- sqrt(colSums(dijv^2))
    cn[[i]] <- sum(dij <= d_max)
  }
  cn
}

#' Calculates number of contacts.
#'
#' @param graph  graph representation of ENM
#' @param d_max cutoff radius to define contacts
#'
#' @return vector of contact numbers
#'
#' @noRd
#'
cn_graph <- function(graph) {
  nsites <- max(c(graph$i, graph$j))
  cn = rep(0,nsites)
  for (i in seq(nsites)) {
    cn[i] <- sum(graph$i == i) + sum(graph$j == i)
  }
  cn
}

#' Reduce 3N x 3N matrix to N x N matrix
#'
#' @noRd
#'
reduce_matrix <- function(full.matrix) {
  N <- nrow(full.matrix) / 3
  reduced.matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      oi <- 3 * (i - 1) + 1
      ooi <- oi + 2
      oj <- 3 * (j - 1) + 1
      ooj <- oj + 2
      Fij <- full.matrix[oi:ooi, oj:ooj]
      reduced.matrix[i, j] <- sum(diag(Fij))
    }
  }
  return(reduced.matrix)
}

#' Distance from each site to the closest active site residue
#'
#' @noRd
#'
dactive.xyz <- function(xyz, site_active) {
  stopifnot(length(xyz) %% 3 == 0)
  nsites <- length(xyz) / 3

  if (anyNA(site_active)) {
    distance <-  rep(NA, nsites)
  } else {
    site <- seq(nsites)
    is_active_site <- site %in% site_active
    xyz <- my_as_xyz(xyz)

    nsite_active <-  sum(is_active_site)
    distance <- rep(NA, nsites)
    for (j in seq(nsites)) {
      d_active_to_j <-  xyz[, j] - xyz[, is_active_site]
      dim(d_active_to_j) <-  c(3, nsite_active)
      d_active_to_j_norm <-  sqrt(colSums(d_active_to_j^2))
      distance[j] <-  min(d_active_to_j_norm)
    }
  }

  distance
}





#' Set of indices of the 3N-long xyz vector that corresponds to a set of sites
#'
#' @param site is a vector of sites
#' @returns a vector of indices such that xyz[xyz_indices_site(site)] are the xyz coorinates of site
#'
#' @noRd
#'
xyz_indices_site <- function(site) {
  s <- rep(site,each=3)
  i <- (s-1)*3+c(1,2,3)
  i
}

#' Fast calculation of quadratic form
#'
#' @noRd
#'
my_quad_form <- function(x,m,y) {
  ret <- crossprod(y,crossprod(m,x))
  as.numeric(ret)
}


#' Matrix trace
#'
#' @noRd
#'
tr <- function(m) sum(diag(m))

#' Log of determinant of matrix
#'
#' @noRd
#'
logdet <- function(m) {
  log(det(m))
}




#' Bhat distance between two distributions
#' (see Fugglebak 2012)
#'
#' @noRd
#'
dbhat <- function(ca, cb, normalize = F) {
  stopifnot(nrow(ca) == nrow(cb))
  if(normalize) {
    ca <- ca/tr(ca)
    cb <- cb/tr(cb)
  }
  dmat = (ca + cb) / 2
  db <- 1/2 * logdet(dmat) -  1/4 * logdet(ca) - 1/4 * logdet(cb)
  db
}

#' RWSIP similarity between two distributions
#' (see Carnevale 2007)
#'
#' @noRd
#'
rwsip <- function(ca, cb, normalize = F) {
  stopifnot(nrow(ca) == nrow(cb))
  if(normalize) {
    ca <- ca/tr(ca)
    cb <- cb/tr(cb)
  }
  eiga <- eigen(ca)
  eigb <- eigen(cb)
  m <- tcrossprod(eiga$values, eigb$values) * crossprod(eiga$vectors, eigb$vectors)^2 / sum(eiga$values * eigb$values)
  rwsip <- sqrt(sum(m))
  rwsip
}

#' Entropy difference between two distribut2ions
#'
#' @noRd
#'
dh <- function(ca, cb, normalize = F) {
  stopifnot(nrow(ca) == nrow(cb))
  if(normalize) {
    ca <- ca/tr(ca)
    cb <- cb/tr(cb)
  }
  eiga <- eigen(ca)
  eigb <- eigen(cb)
  ha <- 1/2 * sum(log(2 * pi * exp(1) * eiga$values))
  hb <- 1/2 * sum(log(2 * pi * exp(1) * eigb$values))
  dh <- hb - ha
  dh
}




