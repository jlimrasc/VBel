\name{VBel-package}
\alias{VBel-package}
\alias{VBel}
\docType{package}
\title{
\packageTitle{VBel}
}
\description{
\packageDescription{VBel}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{VBel}
\packageIndices{VBel}
Define required data i.e.
}
\author{
\packageAuthor{VBel}

Maintainer: \packageMaintainer{VBel}
}
\references{
\url{https://www.tandfonline.com/doi/abs/10.1080/01621459.2023.2169701}
}
\keyword{ package }
\seealso{
\code{\link[=compute_AEL_Rcpp]{compute_AEL_Rcpp()}} for mix of R and C++ computation of AEL

\code{\link[=compute_GVA_Rcpp]{compute_GVA_Rcpp()}} for mix of R and C++ computation of GVA
}
\examples{

#ansGVARcppPure <- compute_GVA_Rcpp(mu, C_0, h, delthh, delth_logpi, z, lam0, 
#rho, elip, a, T, T2, fullCpp = TRUE)
}
