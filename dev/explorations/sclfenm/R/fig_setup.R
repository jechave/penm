library(here)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
theme_set(theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30")))
COL <- c(sclfenm = "#1b6ca8", lfenm = "#d1495b", keepnet = "#66a182",
         neutral = "grey50")
sv <- function(p, name, w = 7, h = 4.2)
  ggsave(here("figures", name), p, width = w, height = h, dpi = 150)
