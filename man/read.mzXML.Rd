% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.mzXML.R
\name{read.mzXML}
\alias{read.mzXML}
\title{Read a standard mzXML}
\usage{
read.mzXML(filename)
}
\arguments{
\item{filename}{A mzXML file as exported by vendor software.}
}
\value{
A generic R object of class mzXML.
}
\description{
\code{read.mzXML} is a copy from the caMassClass package which is no longer actively on CRAN.
}
\details{
The main task of read.mzXML and \code{\link{write.mzXML}} functions is to extract and save scan data of mzXML files. 
In addition attempt is made to keep all other sections of mzXML file as unparsed XML code, so the data can be extracted latter or saved into new mzXML files.
Those unparsed sections are stored as XML text.
}
