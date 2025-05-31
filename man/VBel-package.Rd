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
\author{
\packageAuthor{VBel}

Maintainer: \packageMaintainer{VBel}
}
\references{
References: Yu, W., & Bondell, H. D. (2023), “Variational Bayes for Fast and Accurate Empirical Likelihood Inference”,
Journal of the American Statistical Association, 119(546), 1089–1101. Page 3 \doi{doi:10.1080/01621459.2023.2169701}
}
\keyword{ package }
\seealso{
\code{\link[=compute_AEL]{compute_AEL()}} for C++ computation of AEL

\code{\link[=compute_GVA]{compute_GVA()}} for C++ computation of GVA

\code{\link[=diagnostic_plot]{diagnostic_plot()}} for verifying convergence of computed GVA data

}
\examples{

#ansGVARcppPure <- compute_GVA(mu, C_0, h, delthh, delth_logpi, z, lam0, rho, 
#elip, a, SGD_iters, AEL_iters)
#diagnostic_plot(ansGVARcppPure)
}
