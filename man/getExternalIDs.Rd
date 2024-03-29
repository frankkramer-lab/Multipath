% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reactome.R
\name{getExternalIDs}
\alias{getExternalIDs}
\title{Get External Database IDs of nodes}
\usage{
getExternalIDs(biopax, nodes, database = NULL)
}
\arguments{
\item{biopax}{The biopax object}

\item{nodes}{The list of internal IDs of the nodes}

\item{database}{The name of the database. If missing, the list of all external ids of all databases will be returned.}
}
\value{
A dataframe with the mappings between the internal and external IDs
}
\description{
Get External Database IDs of nodes
}
\examples{
\dontrun{ 
biopax=readBiopax(pi3k.owl)
getExternalIDs(wntBiopax,c("Protein1","Protein2"),"UniProt")
}
}
