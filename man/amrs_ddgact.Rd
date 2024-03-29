% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutscan_amrs_ddgact.R
\name{amrs_ddgact}
\alias{amrs_ddgact}
\title{Calculate ddgact profile using amrs method}
\usage{
amrs_ddgact(wt, pdb_site_active, mut_model, mut_dl_sigma, mut_sd_min)
}
\arguments{
\item{wt}{is the (wild-type) protein to mutate (an object obtained using \code{set_enm})}

\item{pdb_site_active}{is a list of sites (using their resno in the pdb file) that define an "active site"}

\item{mut_model}{is the mutational model, which may be "lfenm" or "sclfenm"}

\item{mut_dl_sigma}{is the standard deviation of a normal distribution from which edge-length perturbations are picked (LFENM model).}

\item{mut_sd_min}{is integer sequence-distance cutoff, only edges with \code{sdij >= mut_sd_min} are mutated}
}
\value{
a tibble containing the ddg profile
}
\description{
`ddgact` is the mutational change in activation free energy calculated using the lfenm model
It is the energy necessary to deform the mutant's active site to "push it back" to its wild-type conformation.
\code{amrs_ddgact()} calculates a profile of site-specific \code{ddgact} values, where the jth value corresponds to averaging over mutations at site j.
}
\details{
For details see \doi{10.7717/peerj.11330}
}
\examples{
\dontrun{
pdb <- bio3d::read.pdb("2acy")
wt <- set_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5, frustrated = FALSE)
pdb_site_active <- c(10,14,21)
responses <- amrs_ddgact(wt, pdb_site_active,  mut_model = "lfenm", mut_dl_sigma = 0.3, mut_sd_min = 1)
}

}
\seealso{
Other mutscan functions: 
\code{\link{admrs}()},
\code{\link{amrs_ddg}()},
\code{\link{amrs}()},
\code{\link{sdmrs}()},
\code{\link{smrs_ddgact}()},
\code{\link{smrs_ddg}()},
\code{\link{smrs}()}
}
\concept{mutscan functions}
