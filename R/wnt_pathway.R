#' Demo function for Wnt Pathway Views
#'
#' @note
#' This function calls "data", which is supposed to be the parsed DrugBank object. Please refer to loadDBXML().
#' @export
#' @importFrom rBiopaxParser readBiopax listPathways
#' @importFrom mully plot3d
wntpathway<-function(){
  downloadPathway("R-HSA-195721")
  wntBiopax=readBiopax("R-HSA-195721.owl")
  pathwayID=getPathwayID(wntBiopax,"R-HSA-195721")
  wntmully=pathway2Mully(wntBiopax,pathwayID)
  plot3d(wntmully,layers=T,vertex.label=NA,edge.width=5,edge.arrow.size=5)
  #First Plot
  plot3d(wntmully,layers=T,vertex.label=NA,edge.width=5)
  #Second Plot
  plot3d(wntmully,layers=T,vertex.label=V(wntmully)$description,edge.width=5)

  view1=pathwayView(wntmully,"View1")
  view1=addStep(view1,action = "remove",element="layer",name="Rna",trans = T)
  suppressWarnings(plot3d(view1$modified))

  view2=pathwayView(wntmully,"View2")
  view2=addStep(view2,action = "remove",element="layer",name="PhysicalEntity",trans=T)
  suppressWarnings(plot3d(view2$modified))
  
  view3=pathwayView(wntmully,"View3")
  view3=addStep(view3,action = "remove",element="layer",name="Complex",trans=T)
  suppressWarnings(plot3d(view3$modified))
  
  view4=undo(view3,1)

  #Add DrugBank Layer
  #Get Proteins on Protein Layer
  proteinLayer=getLayer(wntmully,"Protein")
  proteinMappings=getExternalIDs(wntBiopax,proteinLayer$name,"UniProt")
  proteinIDs=unique(proteinMappings[,2])
  upkbtodb=getUPKBtoDB(proteinIDs)
  drugIDs=unique(upkbtodb$dbid)
  wntmully=addDBLayer(wntmully,data,drugIDs)
  wntmully=removeLayer(wntmully,"drugs",trans = F)
  newWnt=addDBLayer(mully(name="Wnt Pathway",direct = T),dataNew,drugList = drugIDs)
  newWnt=merge(newWnt,wntmully)
  
  #Merge edges
  upkbtodb=read.csv("mappingsWnt.csv")
  dbtoupkb=getDBtoUPKB(data,drugIDs,proteinIDs)
  lnames=intersect(names(upkbtodb),names(dbtoupkb))
  upkbtodb$source="UniProt"
  dbtoupkb$source="DrugBank"
  relations=merge.data.frame(upkbtodb,dbtoupkb[ ,c("dbid","upid","type","dbproteinid","source")],by=lnames,all=T)
  relations$source=paste(relations$'source.x'," ",relations$'source.y')
  relations=relations[,c("dbid","upid","type","dbproteinid","source")]
  relations=unique(relations[order(relations$dbid),])
  names(proteinMappings)=c("intid","upid")
  #Replace internal IDs
  newRelations=merge.data.frame(proteinMappings,newRelations)
  for (i in 1:dim(relations)[1]) {
    startName=V(newWnt)[which(V(newWnt)$name == newRelations$'dbid'[i])]$name
    endName=V(newWnt)[which(V(newWnt)$name == newRelations$'intid'[i])]$name
    attrList=newRelations[i,]
    attr=as.list(attrList)
    names(attr)=names(newRelations)
    newWnt=mully::addEdge(newWnt,startName,endName,attributes = attr[-2])
  }
  wntmully=addGenesLayer(wntmully,wntBiopax)
  wntmully=addDiseaseLayer(wntmully,wntBiopax)
  return(wntmully)
}