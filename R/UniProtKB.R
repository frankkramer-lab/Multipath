#' Get all proteins' entries from UniProt
#'
#' @param up The UniProt.ws Object
#'
#' @return a dataframe containing the Protein's entries with the ID and Name
#' @export
#'
#'
#' @note
#' Should be preceded by UniProt.ws() to get the UniProt.ws Object
#' @examples
#' \dontrun{ 
#' up=UniProt.ws()
#' allProteins=getAllUPKB(up)
#' }
getAllUPKB <- function(up){
  entriesList=UniProt.ws::keys(up,"UniProtKB")
  return(getUPKBInfo(up,entriesList,c("protein_name")))
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
#' To get the list of possible columns, you can call columns(UniProt.ws())
#' @examples
#' \dontrun{ 
#' up <- UniProt.ws()
#' getUPKBInfo(up,c("Q6ZS62","P14384"),c("protein_name","go")) 
#'} 
#' @importFrom UniProt.ws select 
getUPKBInfo <- function(up,proteins,col){
  upcols=UniProt.ws::columns(up)
  boolval=col %in% upcols
  if(FALSE %in% boolval)
    stop("Please check the names of the columns. To find available column names call UniProt.ws::columns(up)")
  dfProt=data.frame(as.list(c("id","UNIPROTKB",col)),stringsAsFactors = FALSE)
  dfProt=dfProt[-1,]
  names(dfProt)=c("id","UNIPROTKB",col)
  skip=F
  for(i in 1:length(proteins)){
    skip=F
    upID=proteins[i]
    #tryCatch to skip errors for non-UniProtKB IDs
    tryCatch({
      row=UniProt.ws::select(up,columns=col,keys=upID,keytype="UniProtKB")
      row=as.list(gsub(';','',row))
      print(paste(i,"Entry Found for",proteins[i],sep = ":")) 
      names(row)=c("id","UNIPROTKB",col)
      dfProt=rbind(dfProt,row)
    },error=function(e){
        skip=T
        warning(paste("The following Protein ID does not exist and will be skipped:",upID))
      },finally={
        if(skip==T)
          next()
        }
      )
  }
  return(dfProt[,-1])
}

#' Get the interactions of given proteins from UniProt
#'
#' @param up The UniProt.ws Object
#' @param proteins The list of proteins of which the interactions should be retrieved
#'
#' @return A dataframe containing the interactions between the given proteins
#' @export
#' 
#' 
#' @note
#' Should be preceded by UniProt.ws() to get the UniProt.ws Object
#' @examples
#' \dontrun{ 
#' up=UniProt.ws()
#' interactions=getUPKBInteractions(up,c("P02747","P07204","P00734"))
#' }
#' @importFrom stringr str_split
getUPKBInteractions<-function(up,proteins){
  allInteractions=getUPKBInfo(up,proteins,c("cc_interaction"))
  interactions=data.frame(V1=is.character(c()),V2=is.character(c()),stringsAsFactors = F)[-1,]
  rows=dim(allInteractions)[1]
  if(rows==0)
    return(interactions)
  for(i in 1:rows){
    listInter=as.list(str_split(allInteractions[i,2],' '))
    listInter=listInter[[1]]
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
#' @param up The UniProt.ws Object
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
#' up=UniProt.ws()
#' g=mully("UniProt")
#' g=addUPKBLayer(g,up,proteinList=c("P02747","P00734","P07204"),col=c("UniProtKB","protein_name"))
#' }
#' @import mully
#' @importFrom svMisc progress
addUPKBLayer<-function(g,up,proteinList,col=c("UniProtKB","protein_name","organism_name")){
  upmully=addLayer(g,"UniProt")
  if(!"UNIPROTKB"%in%col)
    col=append(c("UNIPROTKB"),col)
  proteins=getUPKBInfo(up,proteinList,col)
  interactions=getUPKBInteractions(up,proteins$UNIPROTKB)

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