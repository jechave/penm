library(here)
## FIGURE 7: the response matrix and its two marginals
source(here("R", "sclfenm.R")); source(here("R", "fig_setup.R"))
d <- readRDS(here("data", "response_influence.rds")); R <- d$R; cn <- d$cn
N <- nrow(R)
long <- data.frame(i = rep(1:N, times = N), j = rep(1:N, each = N), v = as.vector(R))
pa <- ggplot(long, aes(j, i, fill = log10(v + 1e-8))) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", name = expression(log[10]~R[ij])) +
  coord_fixed() +
  labs(title = expression(bold("(a)  Response matrix  ")*R[ij]),
       subtitle = "response AT site i (rows) to mutating site j (columns); 2acy chain A",
       x = "mutated site  j", y = "response site  i")
md <- rbind(data.frame(site=1:N, v=rowSums(R), cn=cn,
              what=sprintf("response  (sum over j)   r = %+.2f", cor(rowSums(R),cn))),
            data.frame(site=1:N, v=colSums(R), cn=cn,
              what=sprintf("influence (sum over i)   r = %+.2f", cor(colSums(R),cn))))
pb <- ggplot(md, aes(cn, v, colour = what)) +
  geom_point(size = 1.9, alpha = .85) +
  geom_smooth(method = "lm", formula = y~x, se = FALSE, linewidth = .6) +
  facet_wrap(~what, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("#1b6ca8", "#d1495b"), guide = "none") +
  labs(title = "(b)  The two marginals, against burial",
       subtitle = "buried sites respond less, but influence more",
       x = "contact number", y = NULL)
sv(pa | pb, "fig7_response_matrix.png", 11, 5)
cat(sprintf("response  vs cn: %+.3f\ninfluence vs cn: %+.3f\nsum check: %.4f = %.4f\n",
    cor(rowSums(R),cn), cor(colSums(R),cn), sum(rowSums(R)), sum(colSums(R))))
