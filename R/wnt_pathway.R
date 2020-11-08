#' Demo function for Wnt Pathway Views
#'
#' @param file The link to the Wnt Pathway biopax file 
#'
#' @export
#' @importFrom rBiopaxParser readBiopax listPathways
#' @importFrom mully plot3d
wntpathway<-function(file){
  wntBiopax=readBiopax(file)
  pathwayID=listPathways(wntBiopax)$id[1]
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
}