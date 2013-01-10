\name{peracf}
\alias{peracf}
\title{ Periodic ACF function
}
\description{
Function \code{peracf}, given an input time series and a specified period \code{T}, computes the periodic correlation coefficients such that
\eqn{ \rho(\nu+\tau,\nu)=\rho(\nu,\tau)}, where \eqn{\nu = 1,\ldots, T}
 are seasons and \eqn{\tau} is lag. For each 
possible pair of \eqn{ t} and \eqn{\tau } confidence limits for 
\eqn{\rho(\nu,\tau)} are also computed using Fisher
 transformation. Procedure \code{peracf}
provides also two important tests: \eqn{ \rho(\nu+\tau,\nu) \equiv \rho(\tau)} and \eqn{\rho(\nu+\tau,\nu) \equiv 0}.
}
\usage{
peracf(x, T, tau, missval, datastr,...)
}
\arguments{
  \item{x}{
input time series, at the begining missing values 
in \code{x} will be treat as zeros and periodic mean will be computed, 
then missing values will be replaced by periodic mean. 
}
  \item{T}{
period of PC-T structure.
}
  \item{tau}{
vector of lag values for which estimation is made.
}
  \item{missval}{
notation for missing values (denoted as NaN).
}
 \item{datastr}{
string name of data for printing.
}
\item{...}{other arguments, that are connected with the plots: \code{prttaus, plottaus, cialpha, typeci, typerho, pchci, pchrho, colci, colrho},
where\cr
 \code{prttaus} is a set of lags for which correlation coefficients are printed; it is a subset of \code{tau}, \cr
 \code{plottaus} is a set of lags for plotting the correlation coefficients 
(one plot per lag); it is a subset of \code{tau},\cr
 \code{cialpha} threshold for confidence interval, \cr
 \code{typeci}/ \code{typerho}, \code{pchci}/ \code{pchrho}, \code{colci}/\code{colrho} 
define the type, plot character and colors of confidence intervals/periodic correlation values.\cr
By default these parameters are fixed to \code{prttaus = seq(1,T/2)}, \code{plottaus = seq(1,T/2)}, \code{cialpha = 0.05}, \code{typeci = "b"}, \code{typerho = "b"}, \code{pchci = 10}, \code{pchrho = 15}, \code{colci = "blue"}, \code{colrho = "red"}.
}
}
\details{
Function \code{peracf} uses three separate procedures:\cr
\code{rhoci()} returns the upper and lower bands defining a \eqn{1 - \alpha} confidence interval for the true values of
 \eqn{ \rho(\nu, \tau)},\cr
\code{rho.zero.test()}  tests whether the estimated correlation coefficients are equal to zeros, \eqn{ \rho(\nu+\tau,\nu) \equiv 0}.\cr
\code{rho.equal.test()}  tests whether the estimated correlation coefficients are equal to each other for all seasons in the period,
 \eqn{ \rho(\nu+\tau,\nu) \equiv \rho(\tau)}.\cr

In the test \eqn{ \rho(\nu+\tau,\nu) \equiv \rho(\tau)}, rejection indicates 
that series is properly PC and is not just an amplitude modulated stationary 
sequence. That is, there exist lags for which \eqn{ \rho(\nu+\tau,\nu)}
are properly periodic.\cr
In the test \eqn{\rho(\nu+\tau,\nu) \equiv 0}, the
 rejection for some \eqn{ \tau \neq 0} indicates the sequence is not PC white noise.
}
\value{
tables of values for each specified lag \eqn{\tau}:
\item{rho(t, tau)}{estimated correlation coefficients.}
\item{lower}{lower bands of confidence intervals.}
\item{upper}{upper bands of confidence intervals.}
\item{nsamp}{number of samples used in each estimation.}
}
\references{
Hurd, H. L., Miamee, A. G., (2007), Periodically Correlated Random Sequences: Spectral Theory and Practice, Wiley InterScience.
}
\author{
Harry Hurd 
}
\seealso{
\code{\link{Bcoeff}}, \code{\link{perpacf}}
}
\examples{ 
data(volumes)
dev.set(which=1)
peracf(t(volumes),24,seq(1,12),NaN,'volumes')
}
\keyword{ peracf }

