#' Generate Multipath Graph from General Data
#'
#' @param name The name of the graph to be generated
#' @param proteinList The list of proteins of which the interactions should be retrieved
#' @param data The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#'
#' @return A mully graph with the added data
#' @export
#' @import mully
#' @importFrom igraph V
multipath<-function(name="Multipath",proteinList=NA,data=NA,drugList=NA){
  g=mully(name,direct=T)
  proteinLayer=F
  drugLayer=F

  if(!is.na(drugList) & !is.na(data)){
    message("Multipath: Drug Layer will be added")
    g=addDBLayer(g,data,drugList)
    drugLayer=T
  }
  if(!is.na(proteinList) & !is.na(up)){
    message("Multipath: Protein List will be added")
    g=addUPKBLayer(g,proteinList)
    proteinLayer=T
  }
  #Add Drug-Protein Relations
  if(drugLayer & proteinLayer){
    updbrelations=getUPKBDBRelations(data,proteinList,drugList)
    message("Multipath: Adding DrugBank UniProt Edges")
    for (i in 1:dim(updbrelations)[1]) {
      progress(i,progress.bar = T)
      startName=V(g)[which(V(g)$name == updbrelations$'dbid'[i])]$name
      endName=V(g)[which(V(g)$name == updbrelations$'upid'[i])]$name
      attrList=updbrelations[i,]
      attr=as.list(attrList)
      names(attr)=names(updbrelations)
      g=mully::addEdge(g,startName,endName,attributes = attr[-2])
    }
    message("Multipath: DONE - Adding DrugBank UniProt Edges")
  }
  return(g)
}