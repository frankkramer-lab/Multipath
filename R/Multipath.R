multipath<-function(name="Multipath",up=NA,proteinList=NA,data=NA,drugList=NA){
  g=mully(name)
  proteinLayer=F
  drugLayer=F
  if(!is.na(proteinList) & !is.na(up)){
    message("Multipath: Protein List will be added")
    g=addUPKBLayer(g,up,proteinList)
    proteinLayer=T
  }
  if(!is.na(drugList) & !is.na(data)){
    message("Multipath: Drug Layer will be added")
    g=addDBLayer(g,data,drugList)
    drugLayer=T
  }
  
  #Add Drug-Protein Relations
  if(drugLayer & proteinLayer){
    updbrelations=getUPKBDBRelations(up,data,proteinList,drugList)
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