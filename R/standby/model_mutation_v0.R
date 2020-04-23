dv0_prot <- function(site, prot, dv0_mean = 0, dv0_sd = 1) {
  graph <- prot$enm$graph
  dv0_graph(site, graph, dv0_mean, dv0_sd)
}
dv0_graph <- function(site, graph, dv0_mean = 0, dv0_sd = 1) {
  # calculate dv0 resulting from adding dv0_ij to each contact.
  # assume dv0_ij ~ N(dv0_mean, dv0_sd)
  mutate_edge = graph$i == site | graph$j == site
  graph <- graph[mutate_edge,]
  n_edges <- nrow(graph)
  graph$dv0 = rnorm(n_edges, dv0_mean, dv0_sd)
  sum(graph$dv0)
}

dv0 <- function(cn, mu, sd) {
  stopifnot(length(cn) == 1)
  dv0_contacts = rnorm(cn, mean = mu, sd = sd)
  sum(dv0_contacts)
}

dv_stress_2 <- function(cn, sd_dl = 0.2,
                        n_eff = 1,
                        dg_thr = 2,
                        beta = beta_boltzmann(),
                        k = 4.5) {
  stopifnot(length(cn) == 1)
  ntrials <- 20
  vwt_vector <- rep(NA, ntrials)
  pfix_vector <- rep(NA, ntrials)
  for (trial in seq(ntrials)) {
    s <- rnorm(cn, 0, sd_dl)
    vwt <- sum(.5 * k * s^2)
    dx <- log_fitness_stability(vwt, dg_thr = dg_thr, beta = beta)
    pfix <- p_fix_moran(dx, n_eff)
    pfix_vector[[trial]] <- pfix
    vwt_vector[[trial]] <- vwt
  }
  vwt <- sample(vwt_vector, 1,  prob = pfix_vector)
  s = rnorm(cn, 0, sd_dl)
  vmut <- sum(.5 * k * s^2)
  dv <- vmut - vwt
  dv
}

dv_stress <- function(cn,  sd_g = .2, k = 4.5, sd_dl = 0.2) {
  stopifnot(length(cn) == 1)
  g = rnorm(cn, 0, sd_g)
  dl = rnorm(cn, 0, sd_dl)
  dv1 = sum(.5 * k * dl^2)
  dv2 = -sum(k * g * dl)
  dv1 + dv2
}

dv0_param <- function(cn_mean, cn_sd, mu = 1.3, sd = 1.7) {
  # Taruriki's (2011) distribution of ddg has mu = 1.3 and sd = 1.7
  mu_edge = mu / cn_mean
  sd_edge = sqrt(  (sd ^ 2 - cn_sd ^ 2 *mu_edge ^ 2) / cn_mean )
  lst(mu_edge = mu_edge, sd_edge = sd_edge)
}
