
# Site profile plots ------------------------------------------------------

plot_msf_vs_cn <- function(prot, d_max) {
  tibble(cn = get_cn(prot, d_max), msf = get_msf_site(prot)) %>%
    ggplot(aes(1 / cn, msf)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()
}

plot_bfactor_vs_msf <- function(df, msf, bfactor) {
  tibble(msf = get_msf_site(prot), bfactor = get_bfactor(prot)) %>%
    ggplot(aes(msf, bfactor)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()
}

plot_msf_vs_mlms <- function(prot) {
  tibble(mlms = get_mlms(prot), msf = get_msf_site(prot)) %>%
    ggplot(aes(1 / mlms, msf)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()
}

plot_site_profiles <- function(prot) {

  df <- tibble(
    site = get_site(prot),
    cn = get_cn(prot, prot$enm_param$d_max),
    mlms = get_mlms(prot),
    msf = get_msf_site(prot),
    bfactor = get_bfactor(prot)
  )

  df %>%
    mutate(rcn = 1/cn,
           rmlms = 1 / mlms) %>%
    pivot_longer(
      cols = c(rcn, rmlms, msf, bfactor),
      names_to = "property",
      values_to = "value"
    ) %>%
    mutate(property = factor(property, levels = c("bfactor", "msf", "rcn", "rmlms"),
                             labels = c("bfactor", "msf", "1 / cn", "1 / mlms"))) %>%
    group_by(property) %>%
    mutate(value = jefuns::mMnorm(value)) %>%
    ggplot(aes(y = property, x = site, height = value)) +
    ggridges::geom_ridgeline_gradient() +
    theme_cowplot() +
    theme(axis.title.y = element_blank())
}


# Mode plots --------------------------------------------------------------






#' Plot mode connectivity
#'
plot_Kn <- function(prot) {
  umat2(prot) %>%
    group_by(mode) %>%
    summarise(Kn = Kn(umat2)) %>%
    ggplot(aes(mode, Kn)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()
}

#' Eigenvalues vs. mode
#'
plot_eigenvalue_vs_mode <- function(prot) {
  tibble(mode = get_mode(prot), eigenvalue = get_evalue(prot)) %>%
  ggplot(aes(mode, eigenvalue)) +
    geom_point() +
    theme_cowplot()
}

#' Eigenvalue spectrum
#'
plot_eigenvalue_spectrum <- function(prot) {
  tibble(mode = get_mode(prot), eigenvalue = get_evalue(prot)) %>%
    ggplot(aes(eigenvalue)) +
    geom_histogram() +
    theme_cowplot()
}

#' MSF vs mode
#'
plot_msf_vs_mode <- function(prot) {
  tibble(mode = get_mode(prot), msf = get_msf_mode(prot)) %>%
    ggplot(aes(mode, msf)) +
    geom_point() +
    theme_cowplot()

}

#' MSF spectrum
#'
plot_msf_mode_spectrum <- function(prot) {
  tibble(mode = get_mode(prot), msf = get_msf_mode(prot)) %>%
    ggplot(aes(msf)) +
    geom_histogram() +
    theme_cowplot()
}



# site by mode matrices ---------------------------------------------------



#' Plot umat2 matrix
#'
plot_umat2  <- function(prot) {
  umat2(prot) %>%
    ggplot(aes(x = mode, y = site, fill =sqrt(umat2))) +
    geom_tile() +
    scale_fill_viridis_b() +
    theme_cowplot() +
    NULL
}

#' Plot some normal modes of prot object
#'
plot_modes  <- function(prot, plot_modes = c(1, 2, 3, 4, 5)) {
  umat2(prot) %>%
    filter(mode %in% plot_modes) %>%
    ggplot(aes(y = factor(mode), x = site, height = sqrt(umat2))) +
    geom_ridgeline_gradient() +
    theme_cowplot() +
    NULL

}

#' Plot msf site by mode matrix
#'
plot_msf_site_mode  <- function(prot) {
  msf_site_mode(prot) %>%
    ggplot(aes(x = mode, y = site, fill = sqrt(msf))) +
    geom_tile() +
    scale_fill_viridis_b() +
    theme_cowplot() +
    NULL
}


# Helper functions --------------------------------------------------------



#' Tile plot of a matrix
#'
plot_matrix <- function(m, row_name = "i", col_name = "j", value_name = "mij") {
  df <- matrix_to_tibble(m)
  df %>%
    ggplot(aes(x = j, y = i, fill = mij)) +
    geom_tile() +
    scale_fill_viridis_b() +
    theme_cowplot() +
    xlab(col_name) +
    ylab(row_name) +
    NULL
}







