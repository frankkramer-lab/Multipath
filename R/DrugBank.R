dbmully=mully("DrugBank",direct = T)
dbmully=addLayer(dbmully,c("drugs","targets"))
drugs=all[[1]]
#Add Drugs' nodes
for (i in 1:dim(drugs)[1]) {
  attr=list(dbid=drugs$primary_key[i],type=drugs$type[i])
  dbmully=addNode(dbmully,nodeName = drugs$name[i],layerName = "drugs",attributes = attr)
}

interactions=unique(all[["drug_interactions"]])
#Add Drugs' interactions
for (i in 1:dim(interactions)[1]) {

  startName=V(dbmully)[which(V(dbmully)$dbid == interactions[i,1])]$name
  endName=V(dbmully)[which(V(dbmully)$dbid == interactions$parent_key[i])]$name
  if(!is.null(startName) & !is.null(endName) & length(startName)!=0 & length(endName)!=0 )
    dbmully=addEdge(dbmully,startName,endName,attributes = list(description=interactions$description[i]))
}


targets=all[["drug_targets"]]
#Add Targets' nodes and relations to Drugs
for (i in 1:dim(targets)[1]) {
  targetName=targets$name[i]
  idTarget=getIDNode(dbmully,nameNode = targetName)
  if(is.null(idTarget)){
    attr=list(beid=targets$id[i],organism=targets$organism[i])
    dbmully=addNode(dbmully,nodeName = targets$name[i],layerName = "targets",attributes = attr)
    idTarget=getIDNode(dbmully,nameNode = targetName)
  }
  nameDrug=V(dbmully)[which(V(dbmully)$dbid == targets$parent_key[i])]$name
  if(!is.null(nameDrug)){
    dbmully=addEdge(dbmully,nameDrug,targetName,attributes = list(type="bind"))
  }
}

drugBank=all
for(i in 1:72){
  any(duplicated(all[[i]]$`drugbank-id`))
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rWikiPathways", version = "3.8")


