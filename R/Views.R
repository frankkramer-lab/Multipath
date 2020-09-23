#' Create an empty view
#'
#' @param g The input graph
#' @param name The name of the view
#'
#' @return An empty view
#' @export
#'
#' @examples
#' view=pathwayView(mully("myMully",T),"View1")
pathwayView<-function(g,name){
  id=UUIDgenerate(TRUE)
  op <- options(digits.secs = 6)
  timestamp=Sys.time()
  options(op)
  #Original and last version Graphs
  original=g
  modified=g
  lastStep=0
  steps=data.frame(stepID=is.integer(c()),action=is.character(c()),element=is.character(c()),name=is.character(c()),V1=is.character(c()),V2=is.character(c()),attributesnames=is.list(c()),attributes=is.list(c()),stringsAsFactors=F)[-1,]
  pathwayView=list("id"=id,"name"=name,"original"=original,"modified"=modified,"lastStep"=lastStep,"steps"=steps,"timestamp"=timestamp,"lastmodified"=timestamp)
  class(pathwayView)=c("pathwayView",class(pathwayView))
  return(pathwayView)
}

#' Print function
#'
#' @param pathwayView The input View to be printed
#'
#' @return
#' @export
#' 
print.pathwayView<-function(pathwayView){
  if(missing(pathwayView))
    stop("Invalid Argument")
  cat("View")
  if(!is.na(pathwayView$name)){
    cat(paste(" -- ",green(pathwayView$name),sep=""))
  }
  cat(paste(" -- ",pathwayView$id,sep=""))
  cat(paste(bold("\nCreated: "),blue(pathwayView$timestamp),sep=""))
  cat(paste(" -- ",bold("last modified: "),blue(pathwayView$lastmodified),sep=""))
  cat(paste(bold("\nOriginal Graph: "),green(pathwayView$original$name),sep=""))
  if(dim(pathwayView$steps)[1]==0)
    cat("\nThis View has no steps.")
  else{
    osteps=capture.output(print.data.frame(pathwayView$steps))
    osteps <- paste(osteps, "\n", sep="")
    cat(blue("\nThe View has the following ",bold(dim(pathwayView$steps)[1])," step(s):\n"))
    cat(osteps)
  }
}

is.pathwayView<-function(pathwayView){
  if(!"pathwayView"%in%class(pathwayView))
    return(FALSE)
  if(is.null(pathwayView$original)
     || is.null(pathwayView$id) || is.na(pathwayView$id)
     || is.null(pathwayView$steps))
    return(FALSE)
  return(TRUE)
}

#' Track a modification of a graph
#'
#' @param pathwayView The input view in which the modification should be saved
#' @param action The type of action to be applied. Can either be "add" or "remove
#' @param element The type of the element to be modified. Can either be "node", "edge", or "layer"
#' @param name The name of the element to be modified. This argument is only mandatory for nodes and edges
#' @param layername The layer name. This argument is only mandatory for action "add" and element "node"
#' @param V1 The start node of an edge. This argument is only mandatory for element "edge"
#' @param V2 The end node of an edge. This argument is only mandatory for element "edge"
#' @param attributes The named list of attributes of the element. This argument is required only for action "add". It is optional for both elements "node" and "edge", but mandatory if the edge alread exists 
#' @param multi A boolean whether to select multi-edges or not. This is only mandatory for action "remove" and element "edge". By default set to FALSE, in which case the attributes of the specified edge should be given
#'
#' @return
#' @export
#'
#' @examples
#' g=mully:::demo()
#' view=pathwayView(g,"View1")
#' view=addStep(view,"remove","layer","")
addStep<-function(pathwayView,action,element,name=NA,layername=NA,V1=NA,V2=NA,attributes=NA,multi=F,trans=T){
  if(missing(pathwayView) || !is.pathwayView(pathwayView) || missing(action) || missing(element)
     || !action%in%c("add","remove") || !element%in%c("node","edge","layer") ){
    stop("This step cannot be applied. Please provide a correct view, action and element.")
  }
  if(element=="edge" && (missing(V1) || missing (V2))){
    stop("This step cannot be applied. Please provide the arguments V1 and V2.")
  }
  if(element%in%c("node","layer") && (missing(name) || is.na(name))){
    stop("This step cannot be applied. Please provide the name of the element.")
  }
  if(element=="node" && action=="add" && (missing(layername) || is.na(layername))){
    stop("This step cannot be applied. Please provide the layer name on which the node should be added.")
  }
  
  #tmp variables
  g=pathwayView$modified
  oldg=pathwayView$modified
  stepID=pathwayView$lastStep+1
  steps=data.frame("stepID"=is.integer(c()),"action"=is.character(c()),"element"=is.character(c()),"name"=is.character(c()),"V1"=is.character(c()),"V2"=is.character(c()),"attributesnames"=is.list(c()),"attributes"=is.list(c()),stringsAsFactors = F)[-1,]
  ########## Action Add ############
  #select case add
  if(action=="add"){
    #Call addLayer
    if(element=="layer"){
      g=addLayer(g,name)
      #Update steps in the view
      row=list("stepID"=stepID,"action"="add","element"="layer","name"=name,"V1"=NA,"V2"=NA,"attributes"=NA)
      steps=rbind(steps,row)
    }
    #Call addNode
    if(element=="node"){
      g=addNode(g,name,layername,attributes)
      #Update steps in the view
      att=as.list(getNodeAttributes(g,name))[-1]
      row=list("stepID"=stepID,"action"="add","element"="node","name"=name,"V1"=NA,"V2"=NA,"attributesnames"=paste(names(att),collapse="$$"),"attributes"=paste(att,collapse="$$"))
      steps=rbind(steps,row)
    }
    #Call addEdge
    if(element=="edge"){
      g=addEdge(g,V1,V2,attributes)
      #Update steps in the view
      row=list("stepID"=stepID,"action"="add","element"="edge","name"=NA,"V1"=V1,"V2"=V2,"attributesnames"=paste(names(attributes),collapse="$$"),"attributes"=paste(attributes,collapse="$$"))
      steps=rbind(steps,row)
    }
   
  }
  ########## Action Remove ############
  #select case remove
  if(action=="remove"){
    #Call removeLayer
    if(element=="layer"){
      g=removeLayer(g,name,trans=trans)
      #Update steps in the view
      nodes=getLayer(oldg,name)$name
      ##Add removed edges
      edgesOld=getEdgeAttributes(oldg)
      deletedEdges=rbind(edgesOld[which(edgesOld$V1 %in% nodes),],edgesOld[which(edgesOld$V2 %in% nodes),])
      rows=dim(deletedEdges)[1]
      i=1
      while(i<=rows){
        att=as.list(deletedEdges[i,])[-1][-1]
        row=list("stepID"=stepID,"action"="remove","element"="edge","name"=NA,"V1"=deletedEdges$V1[i],"V2"=deletedEdges$V2[i],"attributesnames"=paste(names(deletedEdges[-1][-1]),collapse="$$"),attributes=paste(att,collapse="$$"))
        steps=rbind(steps,row)
        i=i+1
      }
      
      ##Add removed nodes
      idLayer=oldg$layers$ID[which(oldg$layers$Name==name)]
      nodes=getNodeAttributes(oldg)
      nodes=nodes[which(nodes$n==idLayer),]
      rows=dim(nodes)[1]
      i=1
      while(i<=rows){
        l=as.list(nodes[i,])[-1]
        row=list("stepID"=stepID,"action"="remove","element"="node","name"=nodes[i,1],"V1"=NA,"V2"=NA,"attributesnames"=paste(names(l),collapse="$$"),"attributes"=paste(l,collapse="$$"))
        steps=rbind(steps,row)
        i=i+1
      }
      
      ##Add removed layer
      row=list("stepID"=stepID,"action"="remove","element"="layer","name"=name,"V1"=NA,"V2"=NA,"attributesnames"=NA,"attributes"=NA)
      steps=rbind(steps,row)
    }
    
    #Call removeNode
    if(element=="node"){
      g=removeNode(g,name,trans=trans)
      #Add removed edges
      oldedges=getEdgeAttributes(oldg,name)
      rows=dim(oldedges)[1]
      i=1
      while(i<=rows){
        l=as.list(oldedges[i,])[-1][-1]
        row=list("stepID"=stepID,"action"="remove","element"="edge","name"=NA,"V1"=oldedges$V1[i],"V2"=oldedges$V2[i],"attributesnames"=paste(names(l),collapse="$$"),"attributes"=paste(l,collapse="$$"))
        steps=rbind(steps,row)
        i=i+1
      }
      #Add added transitive edges
      edges=getEdgeAttributes(g)
      if("via"%in%names(edges)){
        transedges=edges[which(edges$type=="trans"),]
        transedges=transedges[which(transedges$via==name),]
        rows=dim(transedges)[1]
        i=1
        while(i<=rows){
          l=as.list(transedges[i,])[-1][-1]
          row=list("stepID"=stepID,"action"="add","element"="edge","name"=NA,"V1"=transedges$V1[i],"V2"=transedges$V2[i],"attributesnames"=paste(names(l),collapse="$$"),"attributes"=paste(l,collapse="$$"))
          steps=rbind(steps,row)
          i=i+1
        }
      }
            #Update steps in the view
      l=as.list(getNodeAttributes(oldg,name))[-1]
      row=list("stepID"=stepID,"action"="remove","element"="node","name"=name,"V1"=NA,"V2"=NA,"attributesnames"=paste(names(l),collapse="$$"),"attributes"=paste(l,collapse="$$"))
      steps=rbind(steps,row)
    }
    
    #Call removeEdge
    if(element=="edge"){
      g=removeEdge(oldg,nodeStart=V1,nodeDest=V2,attributes,multi)
      #Find removed edges
      edgesold=getEdgeAttributes(oldg,V1,V2)
      edgesnew=getEdgeAttributes(g,V1,V2)
      edges=anti_join(edgesold,edgesnew)
      
      rows=dim(edges)[1]
      i=1
      while(i<=rows)
      {
        l=as.list(edges[i,])[-1][-1]
        row=list("stepID"=stepID,"action"="remove","element"="edge","name"=NA,"V1"=V1,"V2"=V2,"attributesnames"=paste(names(l),collapse="$$"),"attributes"=paste(l,collapse="$$"))
        steps=rbind(steps,row)
        i=i+1
      }
    }
    
  }
  
  #Update the steps
  pathwayView$steps=rbind(pathwayView$steps,steps)
  pathwayView$lastStep=stepID
  #Update the modified version in the view
  pathwayView$modified=g
  
  #Change last modification date
  op <- options(digits.secs = 6)
  timestamp=Sys.time()
  options(op)
  pathwayView$lastmodified=timestamp
  class(pathwayView)=c("pathwayView",class(pathwayView))
  return(pathwayView)
}