#' Prot msf and compare with bfactor and 1/cn
#'
plot_msf_site <- function(prot, d_max) {

  df <- get_profiles_site(prot, d_max)

  p_profiles <- df %>%
    mutate(rcn = 1/cn) %>%
    mutate(rwcn = 1/wcn) %>%
    pivot_longer(
      cols = c(rcn, rwcn, msf, bfactor),
      names_to = "property",
      values_to = "value"
    ) %>%
    mutate(property = factor(property, levels = c("bfactor", "msf", "rcn", "rwcn"),
                             labels = c("bfactor", "msf", "1 / cn", "1 / wcn"))) %>%
    group_by(property) %>%
    mutate(value = jefuns::mMnorm(value)) %>%
    ggplot(aes(y = property, x = site, height = value)) +
    ggridges::geom_ridgeline_gradient() +
    theme_cowplot() +
    theme(axis.title.y = element_blank())


  pcn <- df %>%
    ggplot(aes(1 / cn, msf)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()

  pwcn <- df %>%
    ggplot(aes(1 / wcn, msf)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()

  pbfactor <- df %>%
    ggplot(aes(msf, bfactor)) +
    geom_point() +
    geom_smooth() +
    theme_cowplot()




  p2 = plot_grid(pbfactor, pcn, pwcn, ncol = 3)
  plot_grid(p_profiles, p2, ncol = 1)
}


#' Plot some normal modes of prot object
#'
plot_modes  <- function(prot, plot_modes = c(1, 2, 3, 4, 5)) {
  df <- get_msf_site_mode(prot)

  df %>%
    filter(mode %in% plot_modes) %>%
    ggplot(aes(y = factor(mode), x = site, height = sqrt(msf))) +
    geom_ridgeline_gradient() +
    theme_cowplot() +
    NULL

}

#' Plot all normal modes
#'
plot_umat2  <- function(prot) {
  df <- get_umat2(prot)

  df %>%
    # mutate(mode = as.factor(mode)) %>%
    ggplot(aes(x = site, y = mode,  fill = sqrt(msf))) +
    geom_tile() +
    scale_fill_viridis_b() +
    theme_cowplot() +
    NULL

}


get_umat2 <- function(prot) {
  umat <- prot$enm$umat
  mode <- get_mode(prot)
  umat2 <- data.frame(site = prot$site)
  for (n in seq(1:length(mode))) {
    un <- umat[, n]
    un <- my_as_xyz(un)
    un2_column <- data.frame(mode = colSums(un^2))
    names(un2_column) <- paste0("mode_", mode[n])
    umat2 <- cbind(umat2, un2_column)
  }

  df <- umat2 %>%
    pivot_longer(cols = starts_with("mode"),
                 names_to = "mode",
                 names_prefix = "mode_",
                 values_to = "msf") %>%
    mutate(mode = as.integer(mode))  %>%
    arrange(mode, site) %>%
    select(mode, site, msf)
  df
}



