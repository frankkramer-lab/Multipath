#' Add OMIM Layer with its respective Nodes
#'
#' @param g A Mully graph
#' @param biopax A Biopax object
#'
#'
#' @return g A Mully graph including OMIM nodes and possible edges
#' @export
#'
#' @note  should be preceded by calling  romim::set_key('KEY'). The KEY could be requested via omim's official website.
#' @author Mohammad Al Maaz
#' @examples
#' \dontrun{ 
#' biopax=readBiopax("wnt.owl") #This is a mully graph that contains a protein layer
#' pathwayID=listPathways(biopax)$id[1]
#' g=Multipath::pathway2Mully(biopax,pathwayID)
#' g=addGenesLayer(g,biopax)
#' g=addDiseaseLayer(g,biopax)
#' }
addDiseaseLayer <- function(g, biopax, addGenesLayer = FALSE) {
  if (!isLayer(g, "OMIM"))
    g = addLayer(g, "OMIM")
  if (addGenesLayer == TRUE && !isLayer(g, "KEGGGenes"))
    g = addGenesLayer(g, biopax)
  
  if (isLayer(g, "KEGGGenes")) {
    data = getKeggOmimRelation(g, biopax)
    omimIds = data$OMIM
    keggIds = data$KEGGID
    source = data$source
    attributes = data$attributes
    nbnodes = length(getLayer(g,"OMIM")$id)#count of nodes 
    nbedges = 0 #count of edges 
    nbnodesskip = 0
    #Add Genes' Nodes
    message("Multipath: Adding GENE Nodes")
    nodesToAdd = unique(omimIds)
    omiminternalid = c()
    graphnodes = getNodeAttributes(g)
    diseasenodes = dplyr::select(graphnodes, c("database", "id"))
    diseasenodes = diseasenodes[diseasenodes$database %in% c("OMIM"), ]
    for (i in 1:length(nodesToAdd)) {
      if (isFALSE(nodesToAdd[i] %in% diseasenodes$id)) {
        nbnodes = nbnodes + 1
        g = mully::addNode(
          g,
          nodeName = paste0("OMIM", nbnodes),
          layerName = "OMIM",
          attributes =  list(
            "id" = nodesToAdd[i],
            "database" = "OMIM",
            "description" = attributes[i]
          )
        )
        nodeName = paste0("OMIM", nbnodes)
        omiminternalid = c(omiminternalid, nodeName)
        message(
          paste(
            "Multipath: DONE - Disease Node Added:",
            paste0("OMIM", nbnodes),
            "ID=",
            nodesToAdd[i],
            sep = " "
          )
        )
      }
      else{
        message(paste(
          "Skipping: Node already exists",
          "NodeID=",
          nodesToAdd[i],
          sep = " "
        ))
        nbnodesskip = nbnodesskip + 1
        omiminternalid = c(omiminternalid, NA)
        
      }
    }
    
    message("Added Nodes = ", nbnodes)
    message("Skipped Nodes = ", nbnodesskip)
    data = cbind(data, omiminternalid)
    omiminternalid = data$omiminternalid
    
    #Adding OMIM-KEGG Edges
    
    nbedges = 0 #count of edges added
    startName = data$geneInternalId
    endName = omiminternalid
    
    for (i in 1:length(startName)) {
      if (!is.na(endName[i])) {
        g = mully::addEdge(g,
                           nodeStart = startName[i],
                           nodeDest = endName[i],
                           list("source" = source[i]))
        nbedges = nbedges + 1
        message("Multipath: DONE - Disease Edge Added  ",
                paste(endName[i], startName[i], sep = "<>"))
      }
    }
  }
  #Adding OMIM-Proteins Edges
  nbnodes = length(getLayer(g,"OMIM")$id)
  nbedges = 0 #count of edges added
  df_omimfromUPKB = getUPKBRelatedDiseases(g, biopax)
  for (i in 1:length(df_omimfromUPKB$omimid)) {
    nodeName = which(getLayer(g, "OMIM")$id == df_omimfromUPKB$omimid[i])
    if (length(nodeName) == 0) {
      nbnodes = nbnodes + 1
      g = addNode(
        g,
        nodeName = paste0("OMIM", nbnodes),
        layerName = "OMIM",
        attributes = list("id" = df_omimfromUPKB$omimid[i], "database" =
                            "OMIM")
      )
    }
  }
  
  for (i in 1:nrow(df_omimfromUPKB)) {
    skip = F
    omim_internalid = getLayer(g, "OMIM")$name[which(getLayer(g, "OMIM")$id == df_omimfromUPKB$omimid[i])]
    protein_internalid = df_omimfromUPKB$protein_internalid[i]
    tryCatch({
      g = mully::addEdge(g,
                         nodeStart = omim_internalid,
                         nodeDest = protein_internalid,
                         list("source" = "UniProt;OMIM"))
      nbedges = nbedges + 1
      message(
        "Multipath: DONE - Disease Edge Added  ",
        paste(omim_internalid, protein_internalid, sep = "<>")
      )
    }, error = function(e) {
      skip = T
      warning(paste(
        "The following edge already exists and will be skipped:",
        omim_internalid
      ))
    }, finally = {
      if (skip == T)
        next()
    })
  }
  return(g)
}

#' Get proteins that has a reference to OMIM from a mully graph and a biopax
#'
#' @param g        A Mully graph
#' @param biopax   A biopax file
#'
#' @return         Dataframe of Proteins-OMIM relation
#' @export
#' @note  should be preceded by calling  romim::set_key('KEY'). The KEY could be requested via omim's official website.
#' @author Mohammad Al Maaz
#' @examples
#' \dontrun{ 
#' biopax=readBiopax("wnt.owl") 
#' pathwayID=listPathways(biopax)$id[1]
#' g=Multipath::pathway2Mully(biopax,pathwayID) #This is a mully graph that contains a protein layer
#' getUPKBRelatedDiseases(g,biopax)
#' }
getUPKBRelatedDiseases <- function (g, biopax) {
  allextIDs = getExternalIDs(biopax, getLayer(g, "protein")$name)
  extIDs = unique(allextIDs[which(allextIDs$database == "UniProt"),]$extid)
  omimfupkb = getUPKBInfo(extIDs, col = "xref_mim")
  names(omimfupkb) = c("extid", "omimid")
  df_omimfupkb = data.frame(extid = c(),
                            omimid = c())
  
  for (i in 1:nrow(omimfupkb)) {
    omimid = omimfupkb$omimid[i]
    
    # split the omimid into split of six characters
    split = substring(omimid,";")
    for (j in 1:length(split)) {
      new_row = data.frame(extid = omimfupkb$extid[i], omimid = split[j])
      df_omimfupkb = rbind(df_omimfupkb, new_row)
    }
  }
  protein_extid = allextIDs$extid[match(df_omimfupkb$extid, allextIDs$extid)]
  protein_internalid = allextIDs$id[match(df_omimfupkb$extid, allextIDs$extid)]
  upids = data.frame(id =  protein_internalid, extid = protein_extid)
  df_omimfupkb = cbind(df_omimfupkb, upids$id)
  names(df_omimfupkb) = c("uniprotid", "omimid", "protein_internalid")
  return(df_omimfupkb)
}
