## Figures for the klfenm report. Reads data/results.rds.
## Run: Rscript figures/make_figures.R

suppressMessages(library(here))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

res <- readRDS(here("data", "results.rds"))
sv <- function(p, f, w = 9, h = 4.2) {
  ggsave(here("figures", f), p, width = w, height = h, dpi = 150)
  cat("  ", f, "\n")
}
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        strip.background = element_rect(fill = "grey93", colour = NA))
theme_set(th)

## ---- fig 1: contact flips vs sigma -----------------------------------------
d <- res$flip
p1 <- ggplot(d, aes(sigma, frac_any)) +
  geom_line(colour = "#1b6ca8") + geom_point(colour = "#1b6ca8", size = 2) +
  geom_vline(xintercept = 0.3, linetype = 2, colour = "grey50") +
  annotate("text", x = 0.3, y = 0.1, label = "  working sigma", hjust = 0,
           size = 3, colour = "grey40") +
  scale_x_log10() + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "(a)  How often a mutation rewires the network",
       subtitle = "2acy chain A, 200 mutations per sigma, from the wild type",
       x = "mutation size sigma (log scale)", y = "mutations changing the active set")

d2 <- d %>% select(sigma, broken = mean_broken, formed = mean_formed) %>%
  pivot_longer(-sigma)
p2 <- ggplot(d2, aes(sigma, value, colour = name)) +
  geom_line() + geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c(broken = "#d1495b", formed = "#66a182"), name = NULL) +
  labs(title = "(b)  Contacts broken and formed per mutation",
       subtitle = "Roughly balanced at every sigma: the network rewires, it does not simply thin",
       x = "mutation size sigma (log scale)", y = "mean contacts per mutation")
sv(patchwork::wrap_plots(p1, p2, ncol = 2), "fig1_contact_flips.png", 10, 3.8)

## ---- fig 2: downhill moves vs strain ---------------------------------------
d <- res$downhill
p1 <- ggplot(d, aes(strain_e, 100 * frac_neg)) +
  geom_line(colour = "#1b6ca8") + geom_point(size = 2, colour = "#1b6ca8") +
  labs(title = "(a)  Strain opens a downhill channel",
       subtitle = "200 trial mutations from each state; the founder is relaxed and admits none",
       x = "strain energy of the reference state", y = "% of mutations with dV_min < 0")
p2 <- ggplot(d, aes(strain_e)) +
  geom_line(aes(y = mean_dv, colour = "mean"), linewidth = .6) +
  geom_line(aes(y = min_dv, colour = "minimum"), linewidth = .6) +
  geom_point(aes(y = mean_dv, colour = "mean")) +
  geom_point(aes(y = min_dv, colour = "minimum")) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(mean = "#1b6ca8", minimum = "#d1495b"), name = NULL) +
  labs(title = "(b)  The mean stays uphill; the tail crosses zero",
       subtitle = "The mean is flat in strain; the minimum falls below zero",
       x = "strain energy of the reference state", y = "dV_min")
sv(patchwork::wrap_plots(p1, p2, ncol = 2), "fig2_downhill.png", 10, 3.8)

## ---- fig 3: scans, k(l) live vs frozen -------------------------------------
a <- res$scan_on_agg %>% mutate(model = "k(l) live")
b <- res$scan_off_agg %>% mutate(model = "k frozen")
d <- bind_rows(a, b)
p1 <- ggplot(d, aes(site, dv_min, colour = model)) +
  geom_line(alpha = .85) +
  scale_colour_manual(values = c("k(l) live" = "#1b6ca8", "k frozen" = "#d1495b"), name = NULL) +
  labs(title = "(a)  Energy cost of mutating each site",
       subtitle = "Mean over 8 mutations per site, sigma = 0.3", x = "site", y = "dV_min")
p2 <- ggplot(d, aes(cn, dv_min, colour = model)) +
  geom_point(alpha = .7, size = 1.4) +
  scale_colour_manual(values = c("k(l) live" = "#1b6ca8", "k frozen" = "#d1495b"), name = NULL) +
  labs(title = "(b)  Cost against burial", subtitle = "Contact number of the mutated site",
       x = "contact number", y = "dV_min")
p3 <- ggplot(d, aes(site, dts, colour = model)) +
  geom_line(alpha = .85) + geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  scale_colour_manual(values = c("k(l) live" = "#1b6ca8", "k frozen" = "#d1495b"), name = NULL) +
  labs(title = "(c)  Change in conformational entropy",
       subtitle = "Per site; beta = 1 in ANM units, a convention not a temperature",
       x = "site", y = "dTS")
sv(patchwork::wrap_plots(p1, p2, p3, ncol = 3), "fig3_scan.png", 12, 3.6)

## ---- fig 4: trajectories ---------------------------------------------------
d <- bind_rows(lapply(names(res$traj), function(n)
  res$traj[[n]]$record %>% mutate(nu = n)))
p1 <- ggplot(d, aes(step, v_from_founder, colour = nu)) +
  geom_line() +
  labs(title = "(a)  Energy from the founder", subtitle = "Metropolis on the exact dV_min",
       x = "substitutions", y = "V - V_founder")
p2 <- ggplot(d, aes(step, rmsd, colour = nu)) + geom_line() +
  labs(title = "(b)  Structural drift", subtitle = "RMSD from the founder, after superposition",
       x = "substitutions", y = "RMSD")
p3 <- ggplot(d, aes(step, n_active, colour = nu)) + geom_line() +
  labs(title = "(c)  Size of the contact network",
       subtitle = "The network rewires without collapsing", x = "substitutions", y = "active contacts")
sv(patchwork::wrap_plots(p1, p2, p3, ncol = 3), "fig4_trajectories.png", 12, 3.6)

## ---- fig 5: frustration vs rebuild, THE headline ---------------------------
flat <- bind_rows(lapply(res$rebuild, function(seedlist)
  bind_rows(lapply(seedlist, function(z) data.frame(
    seed = z$seed, subs = z$subs, strain = z$strain, strain_e = z$strain_e,
    d_edges = z$d_edges,
    rmsf_max_full  = max(abs(z$full$rmsf$rel)),
    rmsf_cor_full  = z$full$rmsf$cor,
    rmsf_max_cross = max(abs(z$cross$rmsf$rel)),
    rmsf_max_topo  = max(abs(z$topo$rmsf$rel)),
    worst_ov_full  = min(z$full$mode$overlap_best),
    worst_ov_cross = min(z$cross$mode$overlap_best),
    worst_ov_topo  = min(z$topo$mode$overlap_best),
    rmsip_full = z$full$mode$rmsip, rmsip_cross = z$cross$mode$rmsip,
    rmsip_topo = z$topo$mode$rmsip,
    eig_max_full = max(abs(z$full$mode$eigen_rel_diff)) * 100
  )))))

d <- flat %>%
  select(strain_e, cross = rmsf_max_cross, topology = rmsf_max_topo, full = rmsf_max_full) %>%
  pivot_longer(-strain_e, names_to = "cause", values_to = "rmsf_err")
p1 <- ggplot(d, aes(strain_e, rmsf_err, colour = cause)) +
  geom_point(size = 1.8, alpha = .85) + geom_smooth(se = FALSE, method = "loess",
                                                    formula = y ~ x, linewidth = .6) +
  scale_colour_manual(values = c(cross = "#66a182", topology = "#edae49", full = "#d1495b"),
                      name = NULL) +
  labs(title = "(a)  Worst-site RMSF error, decomposed",
       subtitle = "cross = transverse term only; topology = active set only; full = both",
       x = "strain energy of the state", y = "max |RMSF error| over sites (%)")

d2 <- flat %>%
  select(strain_e, cross = worst_ov_cross, topology = worst_ov_topo, full = worst_ov_full) %>%
  pivot_longer(-strain_e, names_to = "cause", values_to = "ov")
p2 <- ggplot(d2, aes(strain_e, ov, colour = cause)) +
  geom_point(size = 1.8, alpha = .85) + geom_smooth(se = FALSE, method = "loess",
                                                    formula = y ~ x, linewidth = .6) +
  scale_colour_manual(values = c(cross = "#66a182", topology = "#edae49", full = "#d1495b"),
                      name = NULL) +
  ylim(0, 1) +
  labs(title = "(b)  Worst individual normal mode",
       subtitle = "Overlap of the least-preserved mode among the 20 softest",
       x = "strain energy of the state", y = "worst |<u_frust, u_rebuilt>|")

d3 <- flat %>% select(strain_e, `worst mode` = worst_ov_full, RMSIP = rmsip_full) %>%
  pivot_longer(-strain_e)
p3 <- ggplot(d3, aes(strain_e, value, colour = name)) +
  geom_point(size = 1.8, alpha = .85) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = .6) +
  scale_colour_manual(values = c(`worst mode` = "#d1495b", RMSIP = "#1b6ca8"), name = NULL) +
  ylim(0, 1) +
  labs(title = "(c)  Why a subspace score is not enough",
       subtitle = "RMSIP stays high while individual modes are lost",
       x = "strain energy of the state", y = "similarity")
sv(patchwork::wrap_plots(p1, p2, p3, ncol = 3), "fig5_frustration_rebuild.png", 12, 3.8)

## ---- fig 6: one state in detail --------------------------------------------
z <- res$rebuild[[1]][[length(res$rebuild[[1]])]]
d <- data.frame(site = seq_along(z$full$rmsf$rel), rel = z$full$rmsf$rel)
p1 <- ggplot(d, aes(site, rel)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_col(fill = "#d1495b", width = .8) +
  labs(title = sprintf("(a)  Per-site RMSF error, one state (%d substitutions)", z$subs),
       subtitle = sprintf("cor = %.3f, yet the worst site is off by %+.1f%%",
                          z$full$rmsf$cor, z$full$rmsf$rel[which.max(abs(z$full$rmsf$rel))]),
       x = "site", y = "RMSF error of the rebuild (%)")
d2 <- data.frame(mode = seq_along(z$full$mode$overlap_best),
                 ov = z$full$mode$overlap_best)
p2 <- ggplot(d2, aes(mode, ov)) +
  geom_col(fill = "#1b6ca8", width = .8) +
  geom_hline(yintercept = z$full$mode$rmsip, linetype = 2, colour = "#d1495b") +
  annotate("text", x = nrow(d2), y = z$full$mode$rmsip, vjust = -.6, hjust = 1,
           label = sprintf("RMSIP = %.3f", z$full$mode$rmsip), size = 3, colour = "#d1495b") +
  ylim(0, 1) +
  labs(title = "(b)  Per-mode overlap, same state",
       subtitle = "The 20 softest modes, matched by maximum overlap",
       x = "mode", y = "|<u_frustrated, u_rebuilt>|")
sv(patchwork::wrap_plots(p1, p2, ncol = 2), "fig6_one_state.png", 10, 3.8)

cat("figures written\n")
