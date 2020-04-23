frustrate_null_boltzmann <-
  function(p_relaxed, beta) {
    # Add a perturbation along each null_mode to lij
    # The coefficient of each null mode is proportional to it's boltzmann prob
    # The equilibrium structure doesn't change
    stopifnot(!is_null(p_relaxed$enm$u_null))
    u_null <- p_relaxed$enm$u_null
    n_null <- ncol(u_null)
    kij <- diag(p_relaxed$enm$graph$kij)
    kmn <- crossprod(u_null, crossprod(kij, u_null))
    eig <- eigen(kmn)
    cnn <- 1/(beta * eig$values)
    an = cnn * rnorm(n_null, mean = 0, sd = 1)
    sij <- sij_an(u_null, an)
    p_frustrated <- p_relaxed
    p_frustrated$enm$graph$lij = p_relaxed$enm$graph$lij + sij
    # recalculate enm, because kmat depends on frustration
    p_frustrated <- update_enm(p_frustrated)
    p_frustrated
  }

frustrate_null_mode <-
  function(p_relaxed, null_mode, s_sd = 1) {
    # it adds perturbation of size s_sd aong mode null_mode to lij
    # the equilibrium conformation remains the same
    stopifnot(!is_null(p_relaxed$enm$u_null))
    u_null <- p_relaxed$enm$u_null
    n_null <- ncol(u_null)
    an = rep(0, n_null)
    an[[null_mode]] <- s_sd
    sij <- sij_an(u_null, an)
    p_frustrated <- p_relaxed
    p_frustrated$enm$graph$lij = p_relaxed$enm$graph$lij + sij
    # recalculate enm, because kmat depends on frustration
    p_frustrated <- update_enm(p_frustrated)
    p_frustrated
  }

frustrate_enm <-
  function(p_relaxed,
           s_sd = .3,
           s_mean = 0,
           sd_min = 2,
           max_iter = 10) {
    # it adds sij to the lij of a relaxed ENM in such a way that
    # the equilibrium conformation remains the same
    stopifnot(!is_null(p_relaxed$enm$u_null))
    u_null <- p_relaxed$enm$u_null
    n_null <- ncol(u_null)


    nmodes = length(p_relaxed$enm$mode)
    for (i in seq(max_iter)) {
      iter = i
      an = rnorm(n_null, mean = s_mean, sd = s_sd)
      sij <- sij_an(u_null, an)
      p_frustrated <- p_relaxed
      p_frustrated$enm$graph$lij = p_relaxed$enm$graph$lij + sij
      # recalculate enm, because kmat depends on frustration
      p_frustrated <- update_enm(p_frustrated)
      nmodes_frustrated <- length(p_frustrated$enm$mode)
      if (nmodes_frustrated == nmodes) {
        break # ok hessian (no negative evalues)
      }
    }

    if (iter == max_iter) print("Warning: saddle point")

    p_frustrated

  }



sij_an <- function(u_null, an) {
  sij <- t(t(u_null) * an)
  sij <- rowSums(sij)
  sij
}


u_null <- function(p_relaxed) {
  # calculate u_null a nodes x edges matrix whose 0-eigenvectors
  # perturb lij of a network without changing its equilibrium structure.
  nsites <- p_relaxed$nsites
  graph <- p_relaxed$enm$graph
  n_edges <- nrow(graph)
  eij <- p_relaxed$enm$eij
  k_node_edge = array(0, dim = c(3, nsites, n_edges))
  for (edge in seq(n_edges)) {
    i <- graph$i[[edge]]
    j <- graph$j[[edge]]
    kij <- graph$kij[[edge]]
    k_node_edge[, i, edge] <- kij * eij[edge, ]
    k_node_edge[, j, edge] <- -kij * eij[edge, ]
  }
  dim(k_node_edge) = c(3 * nsites, n_edges)
  k_edge_edge = crossprod(k_node_edge, k_node_edge)
  eig <- eigen(k_edge_edge)
  null_modes <- near(eig$values, 0)
  u_null <- eig$vectors[, null_modes]
  u_null
}




