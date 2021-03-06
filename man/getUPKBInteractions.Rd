% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UniProtKB.R
\name{getUPKBInteractions}
\alias{getUPKBInteractions}
\title{Get the interactions of given proteins from UniProt}
\usage{
getUPKBInteractions(up, proteins)
}
\arguments{
\item{up}{The UniProt.ws Object}

\item{proteins}{The list of proteins of which the interactions should be retrieved}
}
\value{
A dataframe containing the interactions between the given proteins
}
\description{
Get the interactions of given proteins from UniProt
}
\note{
Should be preceded by UniProt.ws() to get the UniProt.ws Object
}
\examples{
\dontrun{ 
up=UniProt.ws()
interactions=getUPKBInteractions(up,c("P02747","P07204","P00734"))
}
}
