
# Site profile plots ------------------------------------------------------

plot_msf_vs_cn <- function(prot) {
  tibble(cn = get_cn(prot), msf = get_msf_site(prot)) %>%
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
    cn = get_cn(prot),
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
  get_umat2(prot) %>%
    matrix_to_tibble(row_name = "site", col_name = "mode", value_name = "umat2") %>%
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
  get_umat2(prot) %>%
    plot_matrix(row_name = "site", col_name = "mode", value_name = "Uin^2")
}

#' Plot some normal modes of prot object
#'
plot_modes  <- function(prot, plot_modes = c(1, 2, 3, 4, 5)) {
  get_umat2(prot) %>%
    matrix_to_tibble(row_name = "site", col_name = "mode", value_name = "umat2") %>%
    filter(mode %in% plot_modes) %>%
    ggplot(aes(y = factor(mode), x = site, height = sqrt(umat2))) +
    geom_ridgeline_gradient() +
    theme_cowplot() +
    NULL

}

#' Plot msf site by mode matrix
#'
plot_msf_site_mode  <- function(prot) {
  get_msf_site_mode(prot) %>%
    plot_matrix(row_name = "site", col_name = "mode", value_name = "MSF")
}

#' Plot rmsf site by mode matrix
#'
plot_rmsf_site_mode  <- function(prot) {
  get_msf_site_mode(prot) %>%
    sqrt() %>%
    plot_matrix(row_name = "site", col_name = "mode", value_name = "RMSF")
}



# Site by site matrices ---------------------------------------------------


plot_rho_matrix <- function(prot) {
  get_rho_matrix(prot) %>%
    plot_matrix(row_name = "i", col_name = "j", value_name = "rho")
}

plot_reduced_cmat <- function(prot) {
  get_reduced_cmat(prot) %>%
    plot_matrix(row_name = "i", col_name = "j", value_name = "Cij")
}

plot_reduced_kmat <- function(prot) {
  get_reduced_kmat(prot) %>%
    plot_matrix(row_name = "i", col_name = "j", value_name = "Kij")
}



# Helper functions --------------------------------------------------------









