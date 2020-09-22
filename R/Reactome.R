#' Build a mully graph from a given pathway
#'
#' @param biopax The BioPaX object containing the parsed data from an OWL file. This can be obtained using readBiopax(filepath)
#' @param pathwayID The ID of the pathway in the biopax object
#'
#' @return A mully graph built from the given pathway
#' @export
#'
#' @note 
#' This should be preceded by readBiopax(filepath) to obtain the biopax object
#' @examples
#' biopax=readBiopax("pi3k.owl")
#' pi3kmully=pathway2mully(biopax,"pathway1")
pathway2Mully<-function(biopax,pathwayID){
  #Transform created Pathway GraphNEL to igraph object
  pathwaygraph=pathway2Graph(biopax,pathwayID)
  pathwayigraph=graph_from_graphnel(pathwaygraph, name = TRUE, weight = TRUE,unlist.attrs = TRUE)
  
  #Create mully Graph
  properties=getInstanceProperty(biopax,pathwayID)
  pathwayName=properties[1]
  pathwaymully=mully(name=pathwayName,direct = T)  

  #Get Nodes
  listNodes=V(pathwayigraph)$name
  
  #Add Layers
  classes=biopax$dt$class[which(biopax$dt$id%in%listNodes)]
  layers=unique(classes)
  pathwaymully=addLayer(pathwaymully,layers)
  #Add Nodes
  for(i in 1:length(listNodes)){
    nodeAttributes=getInstanceAttributes(biopax,listNodes[i])
    nodeClass=getInstanceClass(biopax,listNodes[i])
    pathwaymully=mully::addNode(pathwaymully,nodeName =  listNodes[i],layerName =  nodeClass,attributes=nodeAttributes[-1])
  } 
  #Add Edges
  edges=as_long_data_frame(pathwayigraph)[,c("from","to")]
  names(edges)=c("V1","V2")
  edgeAttributes=as.data.frame(edge.attributes(pathwayigraph))
  for(i in 1:dim(edges[1])){
    startNode=listNodes[edges$V1[i]]
    endNode=listNodes[edges$V2[i]]
    attributes=edgeAttributes[i,]
    names(attributes)=names(edgeAttributes)
    pathwaymully=mully::addEdge(pathwaymully,startNode,endNode,attributes=attributes)
  }
  return(pathwaymully)
}

getInstanceAttributes<-function(biopax,name){
  XRefs=getXrefAnnotations(biopax,name)
  description=getInstanceProperty(biopax,name)
  l=list("name"=name,id="","description"=description,"database"="")
  for(i in 1:dim(XRefs)[1]){
    string=unlist(str_split(XRefs$annotation[i],":"))
    if(i>1){
      l$id=paste(l$id,",",sep="")
      l$database=paste(l$database,",",sep="")
    }
    l$database=paste(l$database,str_trim(string[1],"both"),sep="")
    l$id=paste(l$id,str_trim(string[2],"both"),sep="")
  }
  return(l)
}