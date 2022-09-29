# Plot comparisons --------------------------------------------------------

#' Plot structural differences along sites
response_site_plot <- function(prot1, prot2) {
  df2 <- delta_structure_df2i(prot1, prot2)
  de2 <- delta_structure_de2i(prot1, prot2)
  "dr2" <- delta_structure_dr2i(prot1, prot2)
  msf <- get_msf_site(prot1)
  df <- tibble(pdb_site = get_pdb_site(prot1), site = get_site(prot1), msf, df2, de2, "dr2")



  plot_site <- df %>%
    pivot_longer(cols = df2:"dr2",
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", ""dr2""),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt("dr2")"))) %>%
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
    pivot_longer(cols = df2:"dr2",
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", ""dr2""),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt("dr2")"))) %>%
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
  # make tibble with data to plot
  df2 <- delta_structure_df2n(prot1, prot2)
  de2 <- delta_structure_de2n(prot1, prot2)
  "dr2" <- delta_structure_dr2n(prot1, prot2)
  msf <- get_msf_mode(prot1)
  df <- tibble( mode = get_mode(prot1), msf, df2, de2, "dr2")

  plot_mode <- df %>%
    pivot_longer(cols = df2:"dr2",
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", ""dr2""),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt("dr2")"))) %>%
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
    pivot_longer(cols = df2:"dr2",
                 names_to = "response_type",
                 values_to = "response") %>%
    group_by(response_type) %>%
    mutate(response = mMnorm(sqrt(response))) %>%
    ungroup() %>%
    mutate(response_type = factor(response_type,
                                  levels = c("df2", "de2", ""dr2""),
                                  labels = c("sqrt(df2)", "sqrt(de2)", "sqrt("dr2")"))) %>%
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


