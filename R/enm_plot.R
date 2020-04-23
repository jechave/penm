#' Prot msf and compare with bfactor and 1/cn
#'
enm_msf_plot <- function(prot, param) {

  df <- enm_profiles(prot, param)

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
enm_modes_plot  <- function(prot, modes = c(1, 2, 3, 4, 5)) {
  plot_modes <- prot$enm$mode %in% modes
  u <- prot$enm$umat[,plot_modes]
  mode <- prot$enm$mode[plot_modes]
  uc2 <- data.frame(site =prot$site)
  for (c in seq(ncol(u))) {
    uc <- u[,c]
    uc <- my_as_xyz(uc)
    uc2_col <- data.frame(mode = colSums(uc^2))
    names(uc2_col) <- paste0("mode_", mode[c])
    uc2 <- cbind(uc2, uc2_col)
  }

  uc2 %>%
    pivot_longer(cols = starts_with("mode"),
                 names_to = "mode",
                 names_prefix = "mode_",
                 values_to = "msf") %>%
    mutate(mode = as.integer(mode)) %>%
    mutate(mode = as.factor(mode)) %>%
    ggplot(aes(y = mode, x = site, height = sqrt(msf))) +
    geom_ridgeline_gradient() +
    theme_cowplot() +
    NULL

}

#' Plot all normal modes
#'
enm_modes_matrix_plot  <- function(prot, modes = prot$enm$mode) {
  plot_modes <- prot$enm$mode %in% modes
  u <- prot$enm$umat[,plot_modes]
  mode <- prot$enm$mode[plot_modes]
  uc2 <- data.frame(site =prot$site)
  for (c in seq(ncol(u))) {
    uc <- u[,c]
    uc <- my_as_xyz(uc)
    uc2_col <- data.frame(mode = colSums(uc^2))
    names(uc2_col) <- paste0("mode_", mode[c])
    uc2 <- cbind(uc2, uc2_col)
  }

  uc2 %>%
    pivot_longer(cols = starts_with("mode"),
                 names_to = "mode",
                 names_prefix = "mode_",
                 values_to = "msf") %>%
    mutate(mode = as.integer(mode)) %>%
    # mutate(mode = as.factor(mode)) %>%
    ggplot(aes(x = site, y = mode,  fill = sqrt(msf))) +
    geom_tile() +
    scale_fill_viridis_b() +
    theme_cowplot() +
    NULL

}


#' Calculate various profiles for prot
#'
enm_profiles <- function(prot, param) {
  site <- prot$site
  pdb_site <- prot$pdb_site
  cn <- cn_xyz(prot$xyz, prot$pdb_site, d_max = param$enm$d_max, sd_min = param$mut$sd_min)
  wcn = wcn_xyz(prot$xyz)
  bfactor <- prot$bfactor
  msf <- msf_prot(prot)

  df <- tibble(site, pdb_site, cn, wcn, bfactor, msf)
  df
}


