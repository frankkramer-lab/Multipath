% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Intraedges.R
\name{getOmimToUPKB}
\alias{getOmimToUPKB}
\title{Get Omim to UPKB relations from OMIM}
\usage{
getOmimToUPKB(omimIds)
}
\arguments{
\item{omimIds}{list of OMIM Ids}
}
\value{
Dataframe of the Proteins-OMIM relation
}
\description{
Get Omim to UPKB relations from OMIM
}
\note{
should be preceded by calling  romim::set_key('KEY'). The KEY could be requested via omim's official website.
}
\examples{
\dontrun{   
library(romim)
romim::set_key('KEY')
OmimToUPKB = getOmimToUPKB(c("611137", "613733", "603816"))
}
}
