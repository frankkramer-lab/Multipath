#' Get Proteins from UniProtKB
#'
#' @param proteins The list of UniProtKB Proteins ID to be retrieved
#' @param col The list of attributes associated to the UniProtKB Entries to be retrieved
#'
#' @return a dataframe containing the protein entries with the selected attributes
#' @export
#'
#'
#' @note
#' To get the list of possible columns, call queryup::returnfields and check the field column"
#' @examples
#' \dontrun{ 
#' getUPKBInfo(c("Q6ZS62","P14384"),c("protein_name","go")) 
#'} 
#' @importFrom queryup get_uniprot_data
getUPKBInfo <- function(proteins,col){
  boolval=col %in% queryup::return_fields$field
  if(FALSE %in% boolval)
    stop(paste(c("The following given:\n",col[which(!boolval)], "\nPlease check the names of the columns. \nTo find available column names, call queryup::returnfields and check the field column")))
  query = list("accession_id" = proteins)
  dfProt = get_uniprot_data(query = query, columns = col)$content
  names(dfProt)=c("entry",col)
  return(dfProt)
}

#' Get the interactions of given proteins from UniProt
#'
#' @param proteins The list of proteins of which the interactions should be retrieved
#'
#' @return A dataframe containing the interactions between the given proteins
#' @export
#' 
#' 
#' @examples
#' \dontrun{ 
#' interactions=getUPKBInteractions(c("P02747","P07204","P00734"))
#' }
#' @importFrom stringr str_split
getUPKBInteractions<-function(proteins){
  allInteractions=getUPKBInfo(proteins,c("cc_interaction"))
  interactions=data.frame(V1=is.character(c()),V2=is.character(c()),stringsAsFactors = F)[-1,]
  rows=dim(allInteractions)[1]
  if(rows==0)
    return(interactions)
  for(i in 1:rows){
    listInter=as.list(str_split(allInteractions[i,2],'; '))
    listInter=listInter[[1]]
    listInter=sapply(listInter, extract_brackets)
    listInter=listInter[which(listInter%in%proteins)]
    if(length(listInter)==0)
      next()
    for(j in 1:length(listInter)){
      entry=c(allInteractions[i,1],listInter[j])
      if(tail(duplicated(rbind(interactions,rev(entry))),1)>0)
        next()
      interactions[dim(interactions)[1]+1,]=entry
    }
  }
  return(interactions)
}

#' Add a protein layer to a mully graph
#'
#' @param g The mully graph
#' @param proteinList The list of UniProt Ids of the proteins to be added
#' @param col The list of attributes associated to the UniProtKB Entries to be retrieved
#'
#' @return The mully graph with the added UniProt layer
#' @export
#'
#' 
#' @note
#' Should be preceded by UniProt.ws() to get the UniProt.ws Object
#' @examples
#' \dontrun{ 
#' g=mully("UniProt")
#' g=addUPKBLayer(g,proteinList=c("P02747","P00734","P07204"),col=c("UniProtKB","protein_name"))
#' }
#' @import mully
#' @importFrom svMisc progress
addUPKBLayer<-function(g,proteinList,col=c("UniProtKB","protein_name","organism_name")){
  upmully=addLayer(g,"UniProt")
  if(!"UNIPROTKB"%in%col)
    col=append(c("UNIPROTKB"),col)
  proteins=getUPKBInfo(proteinList,col)
  interactions=getUPKBInteractions(proteins$UNIPROTKB)
  
  #Add Proteins' Nodes
  message("Multipath: Adding Protein Nodes")
  for (i in 1:dim(proteins)[1]) {
    progress(i,progress.bar = T)
    attrList=proteins[i,]
    attr=as.list(attrList)
    names(attr)=col
    upmully=mully::addNode(upmully,nodeName = proteins$'UNIPROTKB'[i],layerName = "UniProt",attributes = attr[-1])
  }
  message("Multipath: DONE - Protein Nodes Added")
  #Add Proteins' interactions
  message("Multipath: Adding Proteins' Interactions")
  for (i in 1:dim(interactions)[1]) {
    progress(i,progress.bar = T)
    startName=V(upmully)[which(V(upmully)$name == interactions$'V1'[i])]$name
    endName=V(upmully)[which(V(upmully)$name == interactions$'V2'[i])]$name
    if(!is.null(startName) & !is.null(endName) & length(startName)!=0 & length(endName)!=0 )
      upmully=mully::addEdge(upmully,startName,endName,list(source="UniProtKB"))
  }
  message("Multipath: Adding Proteins' Interactions")
  return(upmully)
}

extract_brackets <- function(input) {
  # Check if the input contains square brackets
  if (grepl("\\[.*\\]", input)) {
    # Extract content inside square brackets
    result <- str_extract(input, "(?<=\\[)[^\\]]+(?=\\])")
  } else {
    # Return the input as is
    result <- input
  }
  return(result)
}