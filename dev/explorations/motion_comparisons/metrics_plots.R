# Run from the package root, after metrics_scan.R:
#   Rscript dev/explorations/motion_comparisons/metrics_plots.R
out_dir <- "dev/explorations/motion_comparisons"
scan <- readRDS(file.path(out_dir, "metrics_scan.rds"))
site_mutated <- scan$mutations$site
nsites <- ncol(scan$site$w2_mean)
nmodes <- ncol(scan$mode$w2_mean)
sites <- seq_len(nsites)
modes <- seq_len(nmodes)

scan_mean <- function(level, name) colMeans(scan[[level]][[name]])
unit_sum <- function(x) x / sum(x)
total <- function(level, name) rowSums(scan[[level]][[name]])

colours <- c(mean = "#1f77b4", cov = "#d62728", total = "black",
             kl = "#d62728", jeffreys = "#ff7f0e", bhattacharyya = "#2ca02c", riemann = "#9467bd",
             marginal = "#7f7f7f", w2 = "#1f77b4", frobenius = "#8c564b", dmsf2 = "#e377c2")

plot_parts <- function(level, metric, ylab, main, log = "") {
  x <- if (level == "site") sites else modes
  mean_part <- scan_mean(level, paste0(metric, "_mean"))
  cov_part <- scan_mean(level, paste0(metric, "_cov"))
  both <- mean_part + cov_part
  floor_value <- if (log == "y") min(c(mean_part, cov_part)[c(mean_part, cov_part) > 0]) else 0
  plot(x, both, type = "l", lwd = 2, col = colours["total"], log = log, xlab = level, ylab = ylab, main = main,
       ylim = range(c(both, pmax(mean_part, floor_value), pmax(cov_part, floor_value))))
  lines(x, pmax(mean_part, floor_value), col = colours["mean"])
  lines(x, pmax(cov_part, floor_value), col = colours["cov"])
  legend("topright", c("total", "mean part", "covariance part"), col = colours[c("total", "mean", "cov")],
         lty = 1, lwd = c(2, 1, 1), cex = 0.7, bty = "n")
}

plot_shapes <- function(level, names, labels, cols, main, log = "") {
  x <- if (level == "site") sites else modes
  profiles <- sapply(names, function(name) unit_sum(scan_mean(level, name)))
  matplot(x, profiles, type = "l", lty = 1, col = cols, log = log, xlab = level, ylab = "share of total",
          main = main)
  legend("topright", labels, col = cols, lty = 1, cex = 0.7, bty = "n")
}

plot_joint_vs_marginal <- function(level, joint_names, marginal_name, ylab, main, log = "") {
  x <- if (level == "site") sites else modes
  joint <- Reduce(`+`, lapply(joint_names, function(name) scan_mean(level, name)))
  marginal <- scan_mean(level, marginal_name)
  plot(x, joint, type = "l", lwd = 2, log = log, xlab = level, ylab = ylab, main = main,
       ylim = range(c(joint, marginal)))
  lines(x, marginal, col = colours["marginal"])
  legend("topright", c("joint, decomposed", "sum of marginals"), col = c("black", colours["marginal"]),
         lty = 1, lwd = c(2, 1), cex = 0.7, bty = "n")
}

# --- figure 1 and 2: scan means by site and by mode -------------------------------

for (level in c("site", "mode")) {
  log <- if (level == "mode") "y" else ""
  png(file.path(out_dir, paste0("metrics_by_", level, ".png")), width = 1500, height = 1900, res = 140)
  par(mfrow = c(4, 2), mar = c(4.2, 4.4, 2.2, 1))
  plot_parts(level, "w2", "Å²", "Wasserstein-2 (squared)", log)
  plot_parts(level, "kl", "nats", "KL(mut || wt)", log)
  plot_joint_vs_marginal(level, c("w2_mean", "w2_cov"), "marginal_w2", "Å²", "Wasserstein-2: joint vs marginals", log)
  plot_joint_vs_marginal(level, c("kl_mean", "kl_cov"), "marginal_kl", "nats", "KL: joint vs marginals", log)
  plot_shapes(level, c("kl_cov", "jeffreys_cov", "bhattacharyya_cov", "riemann_cov"),
              c("KL", "Jeffreys", "Bhattacharyya", "Riemannian"),
              colours[c("kl", "jeffreys", "bhattacharyya", "riemann")],
              "covariance part, spectral divergences", log)
  plot_shapes(level, c("w2_cov", "kl_cov", "frobenius_cov", "dmsf2"),
              c("W2", "KL", "Frobenius", "penm dmsf^2"),
              colours[c("w2", "kl", "frobenius", "dmsf2")],
              "covariance part, across families", log)
  plot_shapes(level, c("w2_mean", "kl_mean"), c("W2 (= dr2)", "KL (= beta/2 de2)"),
              colours[c("w2", "kl")], "mean part", log)
  plot(if (level == "site") sites else modes, scan_mean(level, "frobenius_cov"), type = "l", log = log,
       xlab = level, ylab = expression("Å"^4), main = "Frobenius, covariance part", col = colours["frobenius"])
  dev.off()
}

# --- figure 3: three individual mutants, by site ---------------------------------

w2_total <- total("site", "w2_mean") + total("site", "w2_cov")
ranked <- order(w2_total)
examples <- c(low = ranked[round(0.1 * length(ranked))], median = ranked[round(0.5 * length(ranked))],
              high = ranked[length(ranked)])
png(file.path(out_dir, "metrics_examples.png"), width = 1500, height = 1500, res = 140)
par(mfrow = c(3, 2), mar = c(4.2, 4.4, 2.2, 1))
for (label in names(examples)) {
  row <- examples[[label]]
  title_text <- sprintf("%s: site %d, letter %d", label, scan$mutations$site[row], scan$mutations$letter[row])
  for (metric in c("w2", "kl")) {
    mean_part <- scan$site[[paste0(metric, "_mean")]][row, ]
    cov_part <- scan$site[[paste0(metric, "_cov")]][row, ]
    plot(sites, mean_part + cov_part, type = "l", lwd = 2, xlab = "site",
         ylab = if (metric == "w2") "Å²" else "nats",
         main = paste(if (metric == "w2") "W2" else "KL", "-", title_text),
         ylim = range(c(0, mean_part + cov_part)))
    lines(sites, mean_part, col = colours["mean"])
    lines(sites, cov_part, col = colours["cov"])
    abline(v = scan$mutations$site[row], lty = 3)
    legend("topright", c("total", "mean part", "covariance part"), col = colours[c("total", "mean", "cov")],
           lty = 1, lwd = c(2, 1, 1), cex = 0.7, bty = "n")
  }
}
dev.off()

# --- figure 4: global totals across mutants --------------------------------------

kl_mean_total <- total("site", "kl_mean")
kl_cov_total <- total("site", "kl_cov")
w2_mean_total <- total("site", "w2_mean")
w2_cov_total <- total("site", "w2_cov")
kl_marginal_total <- total("site", "marginal_kl")
by_mutated_site <- function(x) tapply(x, site_mutated, mean)

png(file.path(out_dir, "metrics_totals.png"), width = 1500, height = 1500, res = 140)
par(mfrow = c(2, 2), mar = c(4.4, 4.6, 2.2, 1))
plot(w2_mean_total + w2_cov_total, kl_mean_total + kl_cov_total, log = "xy", pch = 16, cex = 0.5,
     xlab = "W2² total (Å²)", ylab = "KL total (nats)", main = "one point per mutant")
plot(kl_mean_total, kl_cov_total, log = "xy", pch = 16, cex = 0.5, col = colours["kl"],
     xlim = range(c(kl_mean_total, w2_mean_total)), ylim = range(c(kl_cov_total, w2_cov_total)),
     xlab = "mean part", ylab = "covariance part", main = "mean vs covariance part, per mutant")
points(w2_mean_total, w2_cov_total, pch = 16, cex = 0.5, col = colours["w2"])
abline(0, 1, lty = 3)
legend("topleft", c("KL (nats)", "W2² (Å²)"), col = colours[c("kl", "w2")], pch = 16, cex = 0.7, bty = "n")
plot(kl_mean_total + kl_cov_total, kl_marginal_total, log = "xy", pch = 16, cex = 0.5,
     xlab = "KL, joint", ylab = "KL, sum of site marginals", main = "joint vs sum of site marginals")
abline(0, 1, lty = 3)
matplot(sites, cbind(unit_sum(by_mutated_site(w2_mean_total + w2_cov_total)),
                     unit_sum(by_mutated_site(kl_mean_total + kl_cov_total)),
                     unit_sum(by_mutated_site(kl_cov_total)),
                     unit_sum(by_mutated_site(w2_cov_total))),
        type = "l", lty = 1, col = c(colours["w2"], colours["kl"], "#ff9896", "#aec7e8"),
        xlab = "mutated site", ylab = "share", main = "total change, by the site mutated")
legend("topright", c("W2 total", "KL total", "KL covariance part", "W2 covariance part"),
       col = c(colours["w2"], colours["kl"], "#ff9896", "#aec7e8"), lty = 1, cex = 0.7, bty = "n")
dev.off()

# --- figure 5: how alike are the scan-mean profiles? -----------------------------

profile_names <- c("w2_mean", "w2_cov", "kl_mean", "kl_cov", "jeffreys_cov", "bhattacharyya_cov", "riemann_cov",
                   "frobenius_cov", "marginal_kl", "marginal_w2", "dmsf2")
spearman_matrix <- function(level) {
  profiles <- sapply(profile_names, function(name) scan_mean(level, name))
  cor(profiles, method = "spearman")
}
png(file.path(out_dir, "metrics_correlations.png"), width = 1800, height = 900, res = 140)
par(mfrow = c(1, 2), mar = c(9, 9, 2.2, 1))
for (level in c("site", "mode")) {
  m <- spearman_matrix(level)
  image(seq_along(profile_names), seq_along(profile_names), m[, rev(seq_along(profile_names))],
        zlim = c(-1, 1), col = hcl.colors(41, "Blue-Red 3"), axes = FALSE, xlab = "", ylab = "",
        main = paste("Spearman between scan-mean profiles, by", level))
  axis(1, seq_along(profile_names), profile_names, las = 2, cex.axis = 0.75)
  axis(2, seq_along(profile_names), rev(profile_names), las = 2, cex.axis = 0.75)
  text(rep(seq_along(profile_names), times = length(profile_names)),
       rep(seq_along(profile_names), each = length(profile_names)),
       sprintf("%.2f", m[, rev(seq_along(profile_names))]), cex = 0.55)
}
dev.off()

# --- numbers for the report -------------------------------------------------------

cat("examples (row, site, letter, W2 total):\n")
for (label in names(examples)) {
  row <- examples[[label]]
  cat(sprintf("  %-6s row %3d  site %2d  letter %d  W2 %.4g\n", label, row, scan$mutations$site[row],
              scan$mutations$letter[row], w2_total[row]))
}
summary_row <- function(x) signif(c(median = median(x), mean = mean(x), min = min(x), max = max(x)), 3)
cat("\nglobal totals over mutants:\n")
print(rbind(
  w2_mean = summary_row(w2_mean_total), w2_cov = summary_row(w2_cov_total),
  kl_mean = summary_row(kl_mean_total), kl_cov = summary_row(kl_cov_total),
  kl_marginal_sites = summary_row(kl_marginal_total)
))
cat("\nfraction of the total in the covariance part, median over mutants: W2",
    signif(median(w2_cov_total / (w2_cov_total + w2_mean_total)), 3),
    " KL", signif(median(kl_cov_total / (kl_cov_total + kl_mean_total)), 3), "\n")
cat("Spearman over mutants, W2 total vs KL total:",
    round(cor(w2_mean_total + w2_cov_total, kl_mean_total + kl_cov_total, method = "spearman"), 3), "\n")
cat("Spearman over mutants, KL mean part vs KL covariance part:",
    round(cor(kl_mean_total, kl_cov_total, method = "spearman"), 3), "\n")
cat("median ratio, KL sum of site marginals / KL joint:",
    signif(median(kl_marginal_total / (kl_mean_total + kl_cov_total)), 3), "\n")
cat("generalised eigenvalues nu over all mutants: min", signif(min(scan$nu_range[, 1]), 3),
    " max", signif(max(scan$nu_range[, 2]), 3),
    "; median per-mutant max |nu - 1|", signif(median(pmax(abs(scan$nu_range[, 1] - 1), abs(scan$nu_range[, 2] - 1))), 3), "\n")
cat("\nSpearman between scan-mean profiles, by site:\n"); print(round(spearman_matrix("site"), 2))
cat("\nSpearman between scan-mean profiles, by mode:\n"); print(round(spearman_matrix("mode"), 2))
