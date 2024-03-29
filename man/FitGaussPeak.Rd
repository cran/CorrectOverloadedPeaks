% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitGaussPeak.R
\name{FitGaussPeak}
\alias{FitGaussPeak}
\title{Extrapolate a flat top peak using a Gauss approach.}
\usage{
FitGaussPeak(
  x,
  y,
  scale_range = c(1, 10),
  steps = 10,
  cutoff = 0.95,
  idx = NULL,
  weight_front = 0.5,
  strip_data = "none",
  account_for_baseline_offset = TRUE,
  method = c("Gauss", "EMG")[1],
  silent = TRUE,
  fix_sd = NULL,
  ...
)
}
\arguments{
\item{x}{A numeric vector, retention times.}

\item{y}{A numeric vector, ion intensities.}

\item{scale_range}{Specifies the expected range for the true peak to exceed the observed, where scale_range=c(1,100) would assume anything between not overloaded and 100-fold overloaded.}

\item{steps}{Specifies a step parameter used to create a sequence within 'scale_range' to test for good fits, higher=more precision, fewer=faster.}

\item{cutoff}{Overloaded peaks will be screwed from Gaussian shape already when approaching detector saturation (DS), cutoff=0.95 ensures that points just before DS will not be used to model fit.}

\item{idx}{If not NULL, 'idx' is expected to specify points to correct explicitly (as a numeric-vector within 1:length(x)).}

\item{weight_front}{A weighting parameter to punish deviations in peak front and tail differently; 0.5=use front/tail equally, 1=use only front, 0=use only tail.}

\item{strip_data}{Use all provided data if 'none' (default). Strip 'front' or 'tail' data in case you observe peak fronting or tailing respectively.}

\item{account_for_baseline_offset}{If TRUE will subtract min(y) from y before fitting parameters.}

\item{method}{The method for peak shape calculation. Can be 'Gauss' or 'EMG' (exponentially modified gauss).}

\item{silent}{For testing purposes some QC-plot will be generated if silent=FALSE.}

\item{fix_sd}{Supply a fix standard deviation (sd) for the peak or leave NULL to estimate sd within function.}

\item{...}{passed to the QC plot function, e.g. 'main' or 'xlab'.}
}
\value{
An annotated plot of the mass spectrum and detailed information within the console (if silent=FALSE) and
the optimal fitted data points (vector of length(y), returned invisible).
}
\description{
\code{FitGaussPeak} will take retention time ('x') and intensity ('y') data and extrapolate all points above
a certain threshold based on further parameters using a Gaussian approach.
}
\details{
This function is mainly used internally (\code{\link{CorrectOverloadedPeaks}}) but can be of value on it's
own to test brute force peak reconstruction given that appropriate base peak chromatograms are available.
}
\examples{
#load test data
data("mzXML_data")
names(mzXML_data)
str(mzXML_data[["scan"]][[1]])
pk <- ModelGaussPeak(height=10^7, width=3, scan_rate=10, e=0, ds=8*10^6, base_line=10^2)
plot(pk, main="Gaussian peak of true intensity 10^7 but cutt off at 8*10^6")
idx <- pk[,"int"]>0.005 * max(pk[,"int"])
tmp <- FitGaussPeak(x=pk[idx,"rt"], y=pk[idx,"int"], silent=FALSE, xlab="RT", ylab="Intensity")

}
