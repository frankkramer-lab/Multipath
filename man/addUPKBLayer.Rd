% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UniProtKB.R
\name{addUPKBLayer}
\alias{addUPKBLayer}
\title{Add a protein layer to a mully graph}
\usage{
addUPKBLayer(
  g,
  proteinList,
  col = c("UniProtKB", "protein_name", "organism_name")
)
}
\arguments{
\item{g}{The mully graph}

\item{proteinList}{The list of UniProt Ids of the proteins to be added}

\item{col}{The list of attributes associated to the UniProtKB Entries to be retrieved}
}
\value{
The mully graph with the added UniProt layer
}
\description{
Add a protein layer to a mully graph
}
\note{
Should be preceded by UniProt.ws() to get the UniProt.ws Object
}
\examples{
\dontrun{ 
g=mully("UniProt")
g=addUPKBLayer(g,proteinList=c("P02747","P00734","P07204"),col=c("UniProtKB","protein_name"))
}
}
