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
#' \dontrun{ 
#' biopax=readBiopax(pi3k.owl)
#' pi3kmully=pathway2mully(biopax,"pathway1")
#' }
#' @import mully
#' @import rBiopaxParser
#' @import igraph
#' @importFrom stringr str_split str_trim
#' @importFrom graph nodeData
pathway2Mully<-function(biopax,pathwayID){
  #Transform created Pathway GraphNEL to igraph object
  pathwaygraph=pathway2Graph(biopax,pathwayID)
  pathwayigraph=graph_from_graphnel(pathwaygraph, name = TRUE, weight = TRUE,unlist.attrs = TRUE)
  
  #Create mully Graph
  properties=getInstanceProperty(biopax,pathwayID)
  pathwayName=properties[1]
  pathwaymully=mully(name=pathwayName,direct = T)  

  #Get Nodes
  ##Empty Graph
  if(length(nodeData(pathwaygraph))==0)
    return(pathwaymully)
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
  for(i in 1:dim(edges)[1]){
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

#' Get External Database IDs of nodes
#'
#' @param biopax The biopax object
#' @param nodes The list of internal IDs of the nodes
#' @param database The name of the database
#'
#' @return A dataframe with the mappings between the internal and external IDs
#' @export
#'
#' @examples
#' \dontrun{ 
#' biopax=readBiopax(pi3k.owl)
#' getExternalIDs(wntBiopax,c("Protein1","Protein2"),"UniProt")
#' }
getExternalIDs<-function(biopax,nodes,database){
  idsRes=data.frame("id"=is.character(c()),"extid"=is.character(c()),stringsAsFactors = F)[-1,]
  for(j in 1:length(nodes)){
    atts=getInstanceAttributes(biopax,nodes[j])
    ids=unlist(str_split(atts$id,","))
    dbs=unlist(str_split(atts$database,","))
    ids=ids[which(dbs==database)]
    idsRes=rbind(idsRes,cbind(rep(nodes[j],length(ids)),ids))
  }
  names(idsRes)=c("id",paste0(database," id"))
  return(idsRes)
}

#' Get internal pathway ID in a BioPAX file
#'
#' @param biopax The biopax object
#' @param reactomeID The Reactome ID of the pathway 
#'
#' @return The internal ID of the pathway in the parsed BioPAX object
#' @note 
#' This should be preceded by readBiopax(filepath) to obtain the biopax object
#' @examples
#' \dontrun{ 
#' biopax=readBiopax(pi3k.owl)
#' id=getPathwayID(biopax,"R-HSA-167057")
#' pi3kmully=pathway2mully(biopax,id)
#' }
#' @importFrom rBiopaxParser listPathways
#' @importFrom stringr str_split
#' @export
getPathwayID<-function(biopax,reactomeID){
  listPathways=listPathways(biopax)
  for(pathway in listPathways$id){
    annots=getInstanceAttributes(biopax,pathway)
    ids=unlist(str_split(annots$id,","))
    dbs=unlist(str_split(annots$database,","))
    id=ids[which(dbs=="Reactome")]
    if(id==reactomeID)
      return(pathway)
  }
  return(NULL)
}

#' Download Reactome Pathways in BioPAX level 2 and 3
#'
#' @param pathwayID The Reactome ID or list of IDs of the pathways to be downloaded. The ID should start with R-HSA-.
#' @param biopaxLevel The BioPAX Level, 2 or 3. By default set to 3.
#' @param destDirectory The Directory in which the Pathway Files should be saved. If missing, the files are saved in the working directory. The Reactome IDs are used to name the files.
#' @param overwrite A Boolean whether to overwrite existing files with the same name.
#'
#' @return The Directory in which the files are saved.
#' @export
#'
#' @examples
#' \dontrun{ 
#' downloadPathway(c("R-HSA-195721","R-HSA-9609507"),biopaxLevel=3,overwrite=T)
#' }
#' @importFrom utils download.file
#' @importFrom RCurl url.exists
#' @importFrom stringr str_split
downloadPathway <-function(pathwayID,biopaxLevel = "3",destDirectory,overwrite = F) {
    if (missing(pathwayID))
      stop("Please provide a Reactome Pathway ID.")
    if (!biopaxLevel %in% c("2", "3"))
      stop("The given Biopax Level is wrong.")
    if (missing(destDirectory)) {
      message("The files will be saved in the Working directory.")
      destDirectory = getwd()
    }
  
    url = paste0("https://reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level",biopaxLevel,"/",sep = "")
    urlcheck="https://reactome.org/content/detail/"
    for (pathway in pathwayID) {
      destFile = paste0(destDirectory, "/", pathway, ".owl", sep = "")
      if (file.exists(destFile)) {
        if (overwrite == T)
          message(paste0("The file ",pathway,".owl already exists and will be overwritten.",sep = ""))
        else{
          message(paste0("The file ",pathway,".owl already exists and will be skipped. If you wish to overwrite it, re-call the function with the argument overwrite set to TRUE.",sep = ""))
          next
        }
          
      }
      
      spl = unlist(str_split(pathway, "R-HSA-"))
      if (length(spl) != 2){
        message(paste0("The following given Reactome ID is wrong and will be skipped: ",pathway,sep=""))
        next
      }
      urlchecktmp=paste0(urlcheck,pathway,sep="")
      urltmp = paste0(url, spl[2], sep = "")
      if (!url.exists(urlchecktmp)){
        message(paste0("The following given Reactome ID is wrong and will be skipped: ",pathway,sep=""))
        next
      }
      download.file(urltmp, destFile)
    }
    message(paste0("The pathways are download and saved in the following directory: ",destDirectory,sep = ""))
    return(destDirectory)
  }