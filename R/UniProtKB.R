data.frame(ID = is.numeric(c()), Name = is.character(c()), NameLower = is.character(c()), stringsAsFactors = FALSE)

up <- UniProt.ws()



getAllUPKB <- function(up){
  entriesList=keys(up,"UNIPROTKB")
  return(getUPKBInfo(up,entriesList,c("PROTEIN-NAMES")))
}

#' Get Proteins from UniProtKB
#'
#' @param up The UniProt.ws Object
#' @param proteins The list of UniProtKB Proteins ID to be retrieved
#' @param col The list of attributes associated to the UniProtKB Entries to be retrieved
#'
#' @return a dataframe containing the protein entries with the selected attributes
#' @export
#'
#'
#' @note
#' Should be preceded by UniProt.ws() to get the UniProt.ws Object
#' @examples
#' getUPKBInfo(c("Q6ZS62","P14384","P40259"),c("PROTEIN-NAMES","DRUGBANK","GO","REACTOME"))
#' To get the list of possible columns, you can call columns(UniProt.ws())
getUPKBInfo <- function(up,proteins,col){
  dfProt=data.frame(as.list(c("UNIPROTKB",col)),stringsAsFactors = FALSE)
  dfProt=dfProt[-1,]
  names(dfProt)=c("UNIPROTKB",col)
  for(i in 1:length(proteins)){
    upID=proteins[i]
    row=select(up,columns=col,keys=upID,keytype="UNIPROTKB")
    dfProt=rbind(dfProt,row)
  }

  return(dfProt)
}


#load Uniprot
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UniProt.ws")
help("UniProt.ws")
##default taxId is Homo Sapiens
up <- UniProt.ws()
entries=keys(up,"UNIPROTKB")
data.frame(ID = is.numeric(c()), Name = is.character(c()), NameLower = is.character(c()), stringsAsFactors = FALSE)

for(i in 1:length(entries)){

}

