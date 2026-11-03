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
    stop("xyz has length ", n3, ", which is not a multiple of 3")
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
#' @returns a vector of indices such that `xyz[xyz_indices_site(site)]` are the xyz coordinates of site
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


#' Stop, or warn, if a protein is not superposed onto a reference
#'
#' Comparisons that use directions (the axes of a site's fluctuations, normal
#' modes) need both proteins in the same frame. The test is how much
#' superposing `prot` onto `reference` would lower their RMSD: for a superposed
#' pair it is zero, up to rounding.
#'
#' @param reference,prot two prots with the same sites
#'
#' @returns nothing; stops if the RMSD would drop by more than 1e-4 Å, warns if
#'   by more than 1e-6 Å
#'
#' @noRd
#'
check_superposed <- function(reference, prot) {
  # RMSD drops, in Å, measured on 2acy: ~1e-16 for superposed pairs (rounding),
  # 3e-7 for a 20-step genm trajectory left unsuperposed, 5e-5 for a mutant
  # rotated by 0.01 degrees, 4e-3 for one rotated by 0.1 degrees
  rmsd_drop_warn <- 1e-6
  rmsd_drop_stop <- 1e-4

  xyz_reference <- get_xyz(reference)
  xyz_prot <- get_xyz(prot)
  all_coordinates <- seq_along(xyz_reference)
  xyz_fitted <- as.vector(bio3d::fit.xyz(fixed = xyz_reference, mobile = xyz_prot,
                                         fixed.inds = all_coordinates,
                                         mobile.inds = all_coordinates))

  # computed here because bio3d::rmsd rounds to 3 decimals, too coarse for this
  difference_as_given <- matrix(xyz_prot - xyz_reference, nrow = 3)   # column k is site k
  difference_fitted <- matrix(xyz_fitted - xyz_reference, nrow = 3)
  rmsd_as_given <- sqrt(mean(colSums(difference_as_given^2)))
  rmsd_fitted <- sqrt(mean(colSums(difference_fitted^2)))
  rmsd_drop <- rmsd_as_given - rmsd_fitted

  if (rmsd_drop > rmsd_drop_stop) {
    stop("the proteins are not superposed: superposing would lower their RMSD by ",
         signif(rmsd_drop, 2), " \u00c5. Superpose them first, with superpose_prot().")
  }
  if (rmsd_drop > rmsd_drop_warn) {
    warning("the proteins are not exactly superposed: superposing would lower their RMSD by ",
            signif(rmsd_drop, 2), " \u00c5. Consider superpose_prot().")
  }
  invisible(NULL)
}


#' Squared overlaps between the modes of a reference and those of a protein
#'
#' Row n is how `reference`'s mode n is distributed over `prot`'s modes. It sums
#' to the share of that mode that lies inside `prot`'s internal space: 1 when the
#' two structures are the same, slightly less when they differ, because each
#' protein's covariance excludes the rotations of its own structure.
#'
#' @param reference,prot two prots with normal modes and the same sites
#'
#' @returns the matrix of squared overlaps, `nmodes x nmodes`; stops if a mode's
#'   share is below 0.9, warns if it is below 0.99
#'
#' @noRd
#'
squared_overlap <- function(reference, prot) {
  # 1 - share, measured on 2acy genm: 7e-7 for a single mutant, 2e-5 for the
  # tip of a 20-step trajectory
  share_warn <- 0.99
  share_stop <- 0.9

  overlap <- crossprod(get_umat(reference), get_umat(prot))
  overlap_squared <- overlap^2
  share <- rowSums(overlap_squared)
  lowest_share <- min(share)

  if (lowest_share < share_stop) {
    stop("a mode of the reference lies only ", signif(100 * lowest_share, 3),
         "% inside the protein's internal space; the proteins are too different to be compared along it")
  }
  if (lowest_share < share_warn) {
    warning("a mode of the reference lies only ", signif(100 * lowest_share, 3),
            "% inside the protein's internal space")
  }
  overlap_squared
}




