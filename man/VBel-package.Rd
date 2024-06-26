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
\code{\link[=compute_AEL]{compute_AEL()}} for choice of R and/or C++ computation of AEL

\code{\link[=compute_GVA]{compute_GVA()}} for choice of R and/or C++ computation of GVA

\code{\link[=diagnostic_plot]{diagnostic_plot()}} for verifying convergence of computed GVA data

}
\examples{

#ansGVARcppPure <- compute_GVA(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, 
#elip, a, iters, iters2, fullCpp = TRUE)
#diagnostic_plot(ansGVARcppPure)
}
