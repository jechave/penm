# Extract the 3x3 covariance block of a single site from a 3N x 3N cmat.
#
# The site-by-site comparisons in R/delta_motion_by_site.R and get_h_site()
# work on these blocks, slicing cmat[, i, , i]. Tests computing the same
# quantities independently use the same input. The full cmat (or
# get_reduced_cmat) is singular by construction, since the six rigid-body modes
# are dropped from the spectrum, so its determinant is zero.
site_block <- function(cmat, i) {
  dim(cmat) <- c(3, nrow(cmat) / 3, 3, nrow(cmat) / 3)
  cmat[, i, , i]
}
