# Plot comparisons --------------------------------------------------------

#' Plot structural differences along sites
response_site_plot <- function(prot1, prot2) {
  df2 <- df2_site(prot1, prot2)
  de2 <- de2_site(prot1, prot2)
  dr2 <- dr2_site(prot1, prot2)
  msf <- msf_site(prot1)
  df <- tibble(pdb_site = prot1$pdb_site, site = prot1$site, msf, df2, de2, dr2)


  plot_site <- df %>%
    pivot_longer(cols = df2:dr2,
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", "dr2"),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt(dr2)"))) %>%
    ggplot(aes(x = site, y = response, color = response_type)) +
    geom_line() +
    ylab("response") +
    scale_color_viridis_d() +
    facet_grid(response_type ~ ., scales = "free_y", switch = "both") +
    # scale_y_log10() +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none") +
    NULL

  plot_msf <- df %>%
    pivot_longer(cols = df2:dr2,
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", "dr2"),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt(dr2)"))) %>%
    ggplot(aes(x =  msf, y = response, color = response_type)) +
    geom_point() +
    geom_smooth() +
    ylab("response") +
    scale_color_viridis_d() +
    facet_grid(response_type ~ ., scales = "free_y", switch = "both") +
    scale_x_log10() +
    scale_y_log10() +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    NULL

  plot_grid(plot_site, plot_msf)
}

#' Plot structural difference along normal modes
response_nm_plot <- function(prot1, prot2) {
  df2 <- df2_nm(prot1, prot2)
  de2 <- de2_nm(prot1, prot2)
  dr2 <- dr2_nm(prot1, prot2)
  msf <- 1 / prot1$enm$evalue
  df <- tibble( mode = prot1$enm$mode, msf, df2, de2, dr2)

  plot_mode <- df %>%
    pivot_longer(cols = df2:dr2,
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", "dr2"),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt(dr2)"))) %>%
    ggplot(aes(x = -mode, y = response, color = response_type)) +
    geom_line() +
    ylab("response") +
    scale_color_viridis_d() +
    facet_grid(response_type ~ ., scales = "free_y", switch = "both") +
    # scale_y_log10() +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none") +
    NULL

  plot_msf <- df %>%
    pivot_longer(cols = df2:dr2,
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", "dr2"),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt(dr2)"))) %>%
    ggplot(aes(x =  msf, y = response, color = response_type)) +
    geom_point() +
    geom_smooth() +
    ylab("response") +
    scale_color_viridis_d() +
    facet_grid(response_type ~ ., scales = "free_y", switch = "both") +
    scale_x_log10() +
    scale_y_log10() +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    NULL

  plot_grid(plot_mode, plot_msf)
}


