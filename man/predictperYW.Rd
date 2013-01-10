\name{predictperYW}
\alias{predictperYW}
\alias{predseries}
\title{
Prediction for PAR model
}
\description{
For known \code{T}, procedure \code{predictperYW} provideds 
the LMS 
forecast of \code{k} future periods (beyond the last point 
in \code{x}) 
of a PAR(p) series; the Yule-Walker method is first use to 
estimate the LMS prediction coefficients using all the 
observed data in \code{x}.\cr
 Additionally, procedure \code{predseries}
 plots the predicted values of the series with real future
 values of the series (provided that such real data is 
available).
}
\usage{
predictperYW(x, T, p, missval, k,...)
predseries(real,x,T,p,...)
}
\arguments{
  \item{x}{
input time series.
}
  \item{T}{
period of PC-T structure.
}
  \item{p}{
order of autoregression, it is assumed constant over time.
}
   \item{missval}{
notation for missing values.
}
  \item{k}{
number of periods for prediction.
}
 \item{real}{
the real future values of \code{x} series (historical data).
}
  \item{...}{other arguments that will be connected with plot: \code{realcol} is a color of konwn values and \code{predcol}
is a color of predicted values on the plot.
By default parameters are fixed to \code{realcol="blue"}, \code{predcol="red"}.
}
}
\value{
function returns and plots two series:
\item{x}{input series with predicted values added.}
\item{new}{only predicted part of the series.}
Function \code{predseries} plots predicted and real values of the series on the same plot.
}
\references{
Box, G. E. P., Jenkins, G. M., Reinsel, G. (1994), Time Series Analysis, 3rd Ed., Prentice-Hall,
Englewood Cliffs, NJ.
\cr
\cr
Brockwell, P. J., Davis, R. A. (1991), Time Series: Theory and Methods, 2nd Ed., Springer: New York.
\cr
\cr
Gladyshev, E. G., (1961), Periodically Correlated Random Sequences, Sov. Math., 2, 385-388.
}
\author{
Wioletta Wojtowicz
}
\examples{
data(volumes)
permest_out<-permest(t(volumes),24, 0.05, NaN,'volumes', pp=0)
xd=permest_out$xd
dev.set(which=1)
predictperYW(xd,24,2,NaN,2)

data(volumes.sep)
dev.set(which=1)
predseries(volumes.sep,t(volumes),24,2)
}
