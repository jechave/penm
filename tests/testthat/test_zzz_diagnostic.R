# TEMPORARY DIAGNOSTIC -- delete once the rwsip NaN on CI is understood.
#
# delta_motion_rwsipi() returns NaN at many sites on the CI macOS runner
# (R 4.6.1, arm64) while returning 1 on the dev machine (R 4.2.1, x86_64).
# It cannot be reproduced locally, so this prints the intermediate values from
# inside rwsip() and fails on purpose, to get them into the CI log.

load(test_path("fixtures", "pdb_2acy_A.rda"))

test_that("DIAGNOSTIC: rwsip intermediates", {
  si <- sessionInfo()
  cat("\n---- environment ----\n")
  cat("R       :", si$R.version$version.string, "\n")
  cat("platform:", si$platform, "\n")
  cat("BLAS    :", si$BLAS, "\n")
  cat("LAPACK  :", si$LAPACK, "\n")

  wt <- set_enm(pdb_2acy_A, node = "ca", model = "ming_wall",
                d_max = 10.5, frustrated = FALSE)
  nsites <- get_nsites(wt)
  cm <- get_cmat(wt)
  dim(cm) <- c(3, nsites, 3, nsites)

  out <- delta_motion_rwsipi(wt, wt)
  cat("\n---- delta_motion_rwsipi(wt, wt) ----\n")
  cat("n NaN :", sum(is.nan(out)), "of", length(out), "\n")
  cat("NaN at:", head(which(is.nan(out)), 20), "\n")
  cat("range of finite values:",
      if (any(is.finite(out))) range(out[is.finite(out)]) else NA, "\n")

  # One compact line per site, so the essentials survive log truncation.
  cat("\nsite sym cplx    minEV        denom        sum(m)       min(m)      rwsip\n")
  for (i in seq_len(nsites)) {
    b <- cm[, i, , i]
    e <- eigen(b)                       # exactly as rwsip calls it
    denom <- sum(e$values * e$values)
    m <- tcrossprod(e$values, e$values) * crossprod(e$vectors, e$vectors)^2 / denom
    r <- suppressWarnings(penm:::rwsip(b, b))
    if (i <= 6L || !isTRUE(all.equal(r, 1))) {
      cat(sprintf("%4d %3s %4s %12.4e %12.4e %12.8f %12.3e %10s\n",
                  i, isSymmetric(b), is.complex(e$values),
                  min(Re(e$values)), denom, sum(Re(m)), min(Re(m)), format(r)))
    }
  }

  expect_true(FALSE)   # deliberate: forces the output into the CI log
})
