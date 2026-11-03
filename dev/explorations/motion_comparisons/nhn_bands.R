# Run from the package root: Rscript dev/explorations/motion_comparisons/nhn_bands.R
out_dir <- "dev/explorations/motion_comparisons"
d <- readRDS(file.path(out_dir, "nhn_scan.rds")); e <- readRDS(file.path(out_dir, "nhn_explain.rds"))
lambda <- d$lambda; n <- length(lambda); nhn_mean <- e$nhn_mean
jump <- which.max(lambda[-1] / lambda[-n])
cat("largest ratio lambda[n+1]/lambda[n] at n =", jump, ":", signif(lambda[jump], 3), "->", signif(lambda[jump + 1], 3), "\n")
cat("modes below the jump:", jump, " above:", n - jump, " (nsites - 1 =", n / 3 + 2 - 1, "backbone bonds)\n")
exp_entropy <- function(w) { p <- w / rowSums(w); pl <- p * log(p); pl[p == 0] <- 0; exp(-rowSums(pl)) }
two_level_nhn <- function(dk, level) { theta <- 0.5 * atan2(2 * abs(dk), abs(outer(level, level, "-"))); w <- tan(theta)^2; diag(w) <- 1; exp_entropy(w) }
model_mean <- function(v, reps = 40, seed = 1) { set.seed(seed); acc <- 0
  for (r in seq_len(reps)) { z <- matrix(rnorm(n * n), n); z <- (z + t(z)) / sqrt(2); acc <- acc + two_level_nhn(z * sqrt(v), lambda) }
  acc / reps }
band <- ifelse(seq_len(n) <= jump, 1, 2)
# coupling sd depends only on the bands of n and m: s_ab, symmetric, three parameters
shape <- function(p) { s <- exp(p); m <- matrix(c(s[1], s[2], s[2], s[3]), 2); m[band, band]^2 }
sse <- function(p) sum((model_mean(shape(p), reps = 15) - nhn_mean)^2)
opt <- optim(log(c(0.04, 0.04, 0.1)), sse, method = "Nelder-Mead", control = list(maxit = 150))
pred <- model_mean(shape(opt$par))
cat("band-constant couplings (3 params) s11, s12, s22 =", signif(exp(opt$par), 3), "\n")
cat("  vs exact: pearson", round(cor(pred, nhn_mean), 3), " spearman", round(cor(pred, nhn_mean, method = "spearman"), 3), "\n")
emp <- d$dk2_mean; diag(emp) <- NA
for (a in 1:2) for (b in a:2) cat("  empirical rms coupling, bands", a, b, ":", signif(sqrt(mean(emp[band == a, band == b], na.rm = TRUE)), 3), "\n")
saveRDS(pred, file.path(out_dir, "nhn_bands_pred.rds"))
png(file.path(out_dir, "nhn_bands.png"), width = 1300, height = 650, res = 140)
par(mar = c(4.2, 4.2, 1, 1))
plot(seq_len(n), nhn_mean, type = "l", lwd = 2, xlab = "mode n", ylab = "mean nhn")
lines(seq_len(n), pred, col = "purple"); lines(seq_len(n), e$s0, col = "red"); abline(v = jump + 0.5, lty = 3)
legend("topright", c("exact", "spectrum + band-constant coupling (3 params)", "spectrum + one constant coupling"),
       col = c("black", "purple", "red"), lty = 1, lwd = c(2, 1, 1), cex = 0.75, bty = "n")
dev.off()
