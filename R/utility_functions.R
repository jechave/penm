
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

cn_xyz <- function(xyz, pdb_site, d_max = 10.5, sd_min = 1) {
  xyz <- my_as_xyz(xyz)
  nsites <- ncol(xyz)
  cn <- rep(NA,nsites)
  for (i in seq(nsites)) {
    sdij <- abs(pdb_site - pdb_site[[i]])
    dijv <- xyz - xyz[,i]
    dij <- sqrt(colSums(dijv^2))
    cn[[i]] <- sum(dij <= d_max & sdij >= sd_min)
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

#' Calculate MSF profile of prot
#'
msf_prot <- function(prot) {
  stopifnot(!is.null(prot$enm$cmat))
  c <- diag(prot$enm$cmat)
  nsites <- length(c) / 3
  dim(c) <- c(3, nsites)
  msf <- colSums(c)
  msf
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


contact_number <- function(d_ij,d_max=12.5) {
  nsites <- ncol(d_ij)
  diag(d_ij) = 0
  cn <- rep(0,nsites)
  for (j in seq(nsites)) {
    cn[j] <- sum(d_ij[,j] <= d_max & d_ij[,j] > 0)
  }
  cn
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



smooth <- function(x, width = 1) {
  # window-smooth function
  if (width == 1) return(x)
  n <- length(x)
  w <- (width - 1)/2
  s <- x
  for (i in seq_along(s)) {
    im <- max(1, i - w)
    ip <- min(n, i + w)
    s[i] = mean(x[im:ip])
  }
  s
}
