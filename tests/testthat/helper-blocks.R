# Extract the 3x3 covariance block of a single site from a 3N x 3N cmat.
#
# This is the shape dbhat/rwsip/dh are actually called with in the package --
# see delta_motion_dbhati/rwsipi/dhi in R/delta_motion_by_site.R, which slice
# cmat[, i, , i] before calling them. Tests must use the same input: the full
# cmat (or get_reduced_cmat) is singular by construction, since the six
# rigid-body modes are dropped from the spectrum, so its determinant is zero
# and logdet has no value on it.
site_block <- function(cmat, i) {
  dim(cmat) <- c(3, nrow(cmat) / 3, 3, nrow(cmat) / 3)
  cmat[, i, , i]
}
