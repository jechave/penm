
beta_boltzmann <- function(R = 1.986e-3, T = 298) 1/(R*T)


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

dr2_xyz <- function(xyz_1, xyz_2) {
  r1 <- my_as_xyz(xyz_1)
  r2 <- my_as_xyz(xyz_2)
  dr = r2 - r1
  colSums(dr^2)
}


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

#' Title Calculate contact number of sites
#'
#' Calculates number of contacts.
#'
#' @param xyz  xyz coordinates of sites (e.g. alpha carbons, or center of mass, etc.)
#' @param d_max cutoff radius to define contacts
#'
#' @return vector of contact numbers
#' @export
#'
#' @examples
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


cn_graph <- function(graph, nsites) {
  cn = rep(0,nsites)
  for (i in seq(nsites)) {
    cn[i] <- sum(graph$i == i) + sum(graph$j == i)
  }
  cn
}




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




distance_to_active <- function(xyz,site_active) {
  stopifnot(length(xyz) %% 3 == 0)
  nsites <- length(xyz)/3

  if (anyNA(site_active)) {
    distance = rep(NA, nsites)
  } else {
    site <- seq(nsites)
    is_active_site <- site %in% site_active
    xyz <- my_as_xyz(xyz)

    nsite.active = sum(is_active_site)
    distance <- rep(NA,nsites)
    for (j in seq(nsites)) {
      d_active_to_j <-  xyz[,j] - xyz[,is_active_site]
      dim(d_active_to_j) <-  c(3,nsite.active)
      d_active_to_j_norm <-  sqrt(colSums(d_active_to_j^2))
      distance[j] <-  min(d_active_to_j_norm)
    }
  }

  distance
}

site_to_ind <- function(site) {
  s <- rep(site,each=3)
  i <- (s-1)*3+c(1,2,3)
  i
}


my_quad_form <- function(x,m,y) {
  ret <- crossprod(y,crossprod(m,x))
  as.numeric(ret)
}
