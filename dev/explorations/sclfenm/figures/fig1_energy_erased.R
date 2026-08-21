library(here)
## FIGURE 1: the refit erases the strain -> dV reported as 0 without the offset
source(here("R", "sclfenm.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5)
set.seed(11); mu <- mutate(wt, 40, sigma = 0.3, model = "sclfenm")

## energy decomposition, as a waterfall
V_stress <- mu$V_stress; V_relax <- mu$V_relax; dV <- mu$dV
d <- data.frame(
  stage = factor(c("1. stress\nat wt geometry","2. relaxation\nrecovered",
                   "3. net cost\n= dV","4. after refit\n(network only)",
                   "5. after refit\n+ saved offset"),
                 levels = c("1. stress\nat wt geometry","2. relaxation\nrecovered",
                            "3. net cost\n= dV","4. after refit\n(network only)",
                            "5. after refit\n+ saved offset")),
  value = c(V_stress, -V_relax, dV, 0, dV),
  kind  = c("term","term","total","WRONG","total"))
p <- ggplot(d, aes(stage, value, fill = kind)) +
  geom_col(width = .62) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_text(aes(label = sprintf("%.3f", value),
                vjust = ifelse(value >= 0, -0.4, 1.3)), size = 3.4) +
  scale_fill_manual(values = c(term = "grey72", total = "#1b6ca8",
                               WRONG = "#d1495b"), guide = "none") +
  labs(title = "The refit erases the strain; the offset puts it back",
       subtitle = "2acy chain A, one mutation at site 40. Bar 4 is what SC-LFENM reports if the energy is not saved before refitting.",
       x = NULL, y = "energy") +
  expand_limits(y = c(-0.5, 1.55))
sv(p, "fig1_energy_erased.png", 8, 4.2)
cat(sprintf("V_stress %.4f  V_relax %.4f  dV %.4f\n", V_stress, V_relax, dV))
cat(sprintf("strain left in refit network: %.2e\n",
    venm(mu$state$r, mu$state$pr, mu$state$k, mu$state$l)))
