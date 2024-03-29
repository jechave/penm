% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enm_analysis.R
\name{get_stress}
\alias{get_stress}
\title{Site-dependent ENM minimum stress energy profile}
\usage{
get_stress(prot)
}
\arguments{
\item{prot}{is a protein object obtained using set_enm()}
}
\value{
a vector of site-dependent stress-energy values
}
\description{
Calculates the sum for each site of the stress energy of each of it's springs at equilibrium
}
\seealso{
Other site profiles: 
\code{\link{get_cn}()},
\code{\link{get_dactive}()},
\code{\link{get_mlms}()},
\code{\link{get_msf_site}()},
\code{\link{get_wcn}()}
}
\concept{site profiles}
