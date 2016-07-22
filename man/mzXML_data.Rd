\name{mzXML_data}
\alias{mzXML_data}
\docType{data}
\title{A test dataset of class \code{mzXML}.}
\description{A test dataset imported by read.mzXML() from a GC-APCI measurement on a Bruker impact II, exported by Compass to an mzXML file.}
\usage{data("mzXML_data")}
\format{
  The full format depends on the source mzXML (different vendors export different fields). However, raw data should be stored within
  list element 'scan' containing sublist elements for 'mass' and 'peaks'. The class-attribute needs to be 'mzXML' such that 
  \code{CorrectOverloadedPeaks} can recognize it.
}
\keyword{datasets}