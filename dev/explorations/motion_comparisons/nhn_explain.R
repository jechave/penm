# Run from the package root: Rscript dev/explorations/motion_comparisons/nhn_explain.R
out_dir <- "dev/explorations/motion_comparisons"
d <- readRDS(file.path(out_dir, "nhn_scan.rds"))
lambda <- d$lambda; n <- length(lambda); idx <- seq_len(n)
nhn_mean <- colMeans(d$nhn)
pt_mean <- colMeans(d$nhn_pt)

exp_entropy <- function(w) { p <- w / rowSums(w); pl <- p * log(p); pl[p == 0] <- 0; exp(-rowSums(pl)) }
two_level_nhn <- function(dk, level) {
  theta <- 0.5 * atan2(2 * abs(dk), abs(outer(level, level, "-")))
  w <- tan(theta)^2; diag(w) <- 1; exp_entropy(w)
}
# mean nhn under random symmetric couplings with variance matrix v (n x n), gaps from lambda
model_mean <- function(v, reps = 40, seed = 1) {
  set.seed(seed)
  acc <- 0
  for (r in seq_len(reps)) {
    z <- matrix(rnorm(n * n), n); z <- (z + t(z)) / sqrt(2)
    acc <- acc + two_level_nhn(z * sqrt(v), lambda)
  }
  acc / reps
}
fit_scale <- function(shape) {
  sse <- function(log_s) sum((model_mean(exp(2 * log_s) * shape) - nhn_mean)^2)
  opt <- optimize(sse, c(-15, 5))
  list(s = exp(opt$minimum), pred = model_mean(exp(2 * opt$minimum) * shape))
}
s0 <- fit_scale(matrix(1, n, n))                       # spectrum only: constant coupling
s2 <- fit_scale(outer(lambda, lambda))                 # spectrum only: coupling ~ s sqrt(lambda_n lambda_m)
s3 <- list(pred = model_mean(d$dk2_mean))              # gaps + empirical mean-square couplings (no fit)

r <- function(a, b) round(c(pearson = cor(a, b), spearman = cor(a, b, method = "spearman")), 3)
cat("mean nhn: min", round(min(nhn_mean), 2), "at n =", which.min(nhn_mean),
    "; max", round(max(nhn_mean), 2), "at n =", which.max(nhn_mean), "\n")
cat("mean nhn at n = 1..5:", round(nhn_mean[1:5], 2), "; last 5:", round(tail(nhn_mean, 5), 2), "\n")
cat("two-level PT with each mutant's own dK vs exact:   ", r(pt_mean, nhn_mean), "\n")
cat("model S3, gaps + empirical mean dK^2 (no fit):     ", r(s3$pred, nhn_mean), "\n")
cat("model S0, spectrum only, constant coupling, s =", signif(s0$s, 3), ":", r(s0$pred, nhn_mean), "\n")
cat("model S2, spectrum only, s sqrt(l_n l_m), s =", signif(s2$s, 3), ":", r(s2$pred, nhn_mean), "\n")

# descriptors of the spectrum
gap_lo <- c(Inf, diff(lambda)); gap_hi <- c(diff(lambda), Inf)
nn_gap <- pmin(gap_lo, gap_hi)
cat("Spearman(mean nhn, -nearest gap):          ", round(cor(nhn_mean, -nn_gap, method = "spearman"), 3), "\n")
cat("Spearman(mean nhn, -nearest relative gap): ", round(cor(nhn_mean, -nn_gap / lambda, method = "spearman"), 3), "\n")
# mode-resolved coupling strength: rms off-diagonal coupling of mode n
dk2 <- d$dk2_mean; diag(dk2) <- NA
coupling_rms <- sqrt(rowMeans(dk2, na.rm = TRUE))
cat("Spearman(mean nhn, rms coupling):          ", round(cor(nhn_mean, coupling_rms, method = "spearman"), 3), "\n")
cat("Spearman(rms coupling, lambda):            ", round(cor(coupling_rms, lambda, method = "spearman"), 3), "\n")
cat("lambda range:", signif(range(lambda), 3), "\n")

saveRDS(list(nhn_mean = nhn_mean, pt_mean = pt_mean, s0 = s0$pred, s2 = s2$pred, s3 = s3$pred,
             nn_gap = nn_gap, coupling_rms = coupling_rms, lambda = lambda),
        file.path(out_dir, "nhn_explain.rds"))

png(file.path(out_dir, "nhn_explain.png"), width = 1500, height = 1100, res = 140)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2, 1))
plot(idx, nhn_mean, type = "l", lwd = 2, xlab = "mode n", ylab = "mean nhn", main = "exact vs models")
lines(idx, pt_mean, col = "grey50", lty = 2); lines(idx, s3$pred, col = "darkgreen")
lines(idx, s0$pred, col = "red"); lines(idx, s2$pred, col = "blue")
legend("topright", c("exact", "PT, own dK", "gaps + mean dK^2", "spectrum, const coupling", "spectrum, s sqrt(l l)"),
       col = c("black", "grey50", "darkgreen", "red", "blue"), lty = c(1, 2, 1, 1, 1), cex = 0.7, bty = "n")
plot(idx, lambda, type = "l", log = "y", xlab = "mode n", ylab = "lambda (log)", main = "spectrum")
plot(idx, nn_gap, type = "l", log = "y", xlab = "mode n", ylab = "nearest gap (log)", main = "nearest-neighbour gap")
lines(idx, nn_gap / lambda * median(lambda), col = "blue")
legend("topleft", c("absolute", "relative x median(lambda)"), col = c("black", "blue"), lty = 1, cex = 0.7, bty = "n")
plot(idx, coupling_rms, type = "l", log = "y", xlab = "mode n", ylab = "rms coupling (log)", main = "mean-square coupling of mode n")
dev.off()
