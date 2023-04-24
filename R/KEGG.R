# Get A Gene using KeggGet
#' Get Kegg gene
#'
#' @param geneList  List of Genes as formatted in KEGG Genes
#'
#' @return A data frame of kegg data
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' genes=  getKeggGene(geneList)
getKeggGene <- function(geneList) {
  if (missing(geneList) || is.null(geneList)) {
    stop("Invalid Arguments")
  }
  df_keggGetInput = data.frame(
    entry = is.character(c()),
    organism = is.character(c()),
    name = is.character(c()),
    idLinks = is.character(c()),
    stringsAsFactors = FALSE
  )[-1, ]
  h = 1
  i = 1
  while (i <= length(geneList)) {
    #Set of 10 elements
    indexEnd = h * 10
    #Check if out of bound
    if (indexEnd > length(geneList))
      indexEnd = length(geneList)
    list_Get =  c(geneList[i:indexEnd])
    list_keggGet = KEGGREST::keggGet(c(geneList[i:indexEnd]))
    genes = transformKeggData(list_keggGet, list_Get)
    df_keggGetInput = merge.data.frame(
      x = df_keggGetInput,
      y = genes,
      all.x = TRUE,
      all.y = TRUE
    )
    i = indexEnd + 1
    h = h + 1
  }
  return(df_keggGetInput)
}

#' Transform data retrieved from KEGG
#'
#' @param list_keggGet list of gene input
#' @param list_Get  list of gene data retrieved from Kegg genes database
#'
#'@note must be precedeed by getkegggenes(geneInput)
#' @return a dataframe of Kegg data
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' genes=  getKeggGene(geneList)
transformKeggData <- function(list_keggGet, list_Get) {
  df_keggGet = data.frame(
    entry = is.character(c()),
    organism = is.character(c()),
    name = is.character(c()),
    idLinks = is.character(c()),
    stringsAsFactors = FALSE
  )[-1,]
  for (m in 1:length(list_keggGet)) {
    #entry = toString(list_keggGet[[m]][["ENTRY"]]) # this returns the entry id in KEGG database.
    entry = toString(list_Get[m])  #This returns the name same as the user entered in the gene list. Better for understanding
    name = toString(list_keggGet[[m]][["NAME"]])
    organism = toString(list_keggGet[[m]][["ORGANISM"]])
    idLinks = toString(unlist(str_split(list_keggGet[[m]][["DBLINKS"]], ": ")))
    df_keggGettemp = data.frame(
      entry = entry,
      organism = organism,
      name = name,
      idLinks = idLinks
    )
    df_keggGet =  merge.data.frame(
      x = df_keggGet,
      y = df_keggGettemp,
      all.x = TRUE,
      all.y = TRUE
    )
  }
  return(df_keggGet)
}

#Get interactions between a kegg gene input and all the databases - returns NA when there is no cross reference to other databases
#' Interaction between KEGG genes and other databases
#'
#' @param keggInput data frame of all data retrieved from kegg genes database
#'
#' @return A data frame of each gene vs interactions with other databases
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' genes = getKeggGene(geneList)
#' genesinteraction = interactionKegg(genes)
interactionKegg <- function(keggInput) {
  df_keggLinks = data.frame(
    entry = is.character(c()),
    organism = is.character(c()),
    stringsAsFactors = FALSE
  )[-1,]
  
  for (i in 1:nrow(keggInput)) {
    entry = keggInput[i, 1]
    organism = keggInput[i, 2]
    idList = keggInput[i, 4]
    links =   (unlist(str_split(idList, ", ")))
    #Get odd located elements
    Names = c("entry", "organism", links[c(TRUE, FALSE)])
    #Get even located elements
    ids = c(entry, organism, links[c(FALSE, TRUE)])
    ids = str_replace(ids, " ", ";")
    #Create named list
    names(ids) = Names
    #Transform to data frame and switch columns/rows
    df_keggGettemp = t(BiocGenerics:::as.data.frame(ids))
    df_keggLinks = merge.data.frame(
      x = df_keggLinks,
      y = df_keggGettemp,
      all.x = TRUE,
      all.y = TRUE
    )
  }
  
  return(df_keggLinks)
}

#simplifies the interaction between a gene input and other databases and shows results in 4 columns (c1 = entry number in KEGG c2= name in other databse c3= other database's name c4= organism name as attribtute  )
#' Simplify the dataframe's structure
#'
#' @param genes A data frame of all genes and their relation to other cross references
#'
#' @note  must be preceeded by interactionKegg(keggInput)
#' @return A data frame of gene and its presence in other databases
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' genes = getKeggGene(geneList)
#' genesinteraction = interactionKegg(genes)
#' genesinteractionsimplified =simplifyInteractionKegg(genesinteraction)
simplifyInteractionKegg <- function(genes) {
  df_simplifyInteraction = data.frame(
    keggEntry = is.character(c()),
    dbId = is.character(c()),
    dbName = is.character(c()),
    attributes = is.character(c()),
    source = is.character(c()),
    stringsAsFactors = FALSE
  )[-1,]
  #  if ("NCBI-GeneID" %in% colnames(genes))
  #    genes = genes[, -which(colnames(genes) %in% c("NCBI-GeneID"))]
  entry = genes[1]
  attributes = genes[2]
  dblinks = genes[3:ncol(genes)]
  for (i in 1:dim(entry)[1]) {
    for (j in 1:ncol(dblinks)) {
      if (!is.na(dblinks[i, j])) {
        idlink =   (unlist(str_split(dblinks[i, j], ";")))
        for (l in 1:length(idlink)) {
          keggEntry  = genes$entry[i]
          dbId =  idlink[l]
          dbName = colnames(dblinks[j])
          attributes = genes$organism[i]
          df_simplifyInteractiontemp = data.frame(
            keggEntry = keggEntry,
            dbId = dbId,
            dbName = dbName,
            attributes = attributes,
            source = "KEGGGENES",
            stringsAsFactors = FALSE
          )
          df_simplifyInteraction =  merge.data.frame(
            x = df_simplifyInteraction,
            y = df_simplifyInteractiontemp,
            all.x = TRUE,
            all.y = TRUE
          )
        }
      }
    }
  }
  return(df_simplifyInteraction)
}

#' Add KEGG Gene Layer with its respective Nodes
#'
#' @param g  An existing Mully graph
#' @param biopax A biopax file
#'
#' @return A Mully graph with the all layers,nodes and edges
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' biopax=readBiopax("wnt.owl") #This is a mully graph that contains a protein layer
#' pathwayID=listPathways(biopax)$id[1]
#' g=Multipath::pathway2Mully(biopax,pathwayID)
#' g=addGenesLayer(g,biopax)
addGenesLayer <- function(g, biopax) {
  g = addLayer(g, "KEGGGenes")
  data = getKeggUpkbRelations(g, biopax)
  genesList = data$keggid
  proteinsList = data$upid
  sourcekegg = data$source.y
  sourcedb = data$source.x
  internalid = data$internalid
  geneattributes = data$keggattributes
  kegginternalid = c()
  nbnodesskip = 0 #count of nodes skipped
  nbnodes = 0 #count of nodes added
  nbedges = 0 #count of edges added
  #Add Genes' Nodes
  message("Multipath: Adding GENE Nodes")
  nodesToAdd = genesList
  for (i in 1:length(nodesToAdd)) {
    graphnodes = getNodeAttributes(g)
    genenodes = dplyr::select(graphnodes, c("database", "id"))
    genenodes = genenodes[genenodes$database %in% c("KEGGGenes"), ]
    if (isFALSE(nodesToAdd[i] %in% genenodes$id)) {
      nbnodes = nbnodes + 1
      g = mully::addNode(
        g,
        nodeName = paste0("KEGGGene", nbnodes),
        layerName = "KEGGGenes",
        attributes =  list(
          "id" = nodesToAdd[i],
          "database" = paste(sourcekegg[i], sourcedb[i], sep = ", "),
          "description" = geneattributes[i]
        )
      )
      nodeName = paste0("KEGGGene", nbnodes)
      kegginternalid = c(kegginternalid, nodeName)
      message(paste(
        "Multipath: DONE - GENE Node Added:",
        paste0("KEGGGene", nbnodes),
        "ID=",
        nodesToAdd[i],
        sep = " "
      ))
    }
    else{
      message(paste(
        "Skipping: Node already exists",
        "NodeID=",
        nodesToAdd[i],
        sep = " "
      ))
      nbnodesskip = nbnodesskip + 1
      kegginternalid = c(kegginternalid, NA)
      
    }
  }
  message("Added Nodes = ", nbnodes)
  message("Skipped Nodes = ", nbnodesskip)
  data = cbind(data, kegginternalid)
  kegginternalid = data$kegginternalid
  
  #Adding Edges
  nbedges = 0 #count of edges added
  startName = data$internalid
  endName = data$kegginternalid
  for (i in 1:length(genesList)) {
    if (!is.na(endName[i])) {
      g = mully::addEdge(g,
                         nodeStart = startName[i],
                         nodeDest = endName[i],
                         list("source" = paste(sourcekegg[i], sourcedb[i], sep =
                                                 ",")))
      nbedges = nbedges + 1
      message("Multipath: DONE - GENE Edge Added  ",
              paste(endName[i], startName[i], sep = "<>"))
    }
    
  }
  return(g)
}

#' Get proteins that has a reference to KEGG Genes from a mully graph and a biopax
#'
#' @param g        A Mully graph
#' @param biopax   A biopax file
#'
#' @return         List of Protein
#' @export
#'
#' @examples
#'
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' relatedGenes = getRelatedGenes(g,biopax)
getRelatedGenes <- function (g, biopax) {
  allextIDs = getExternalIDs(biopax, getLayer(g, "protein")$name)
  extIDs = unique(allextIDs[which(allextIDs$database == "UniProt"),]$extid)
  #genesfupkb = getUPKBInfo(UniProt.ws(), extIDs, col = c("UNIPROTKB", "KEGG"))
  genesfupkb = getUPKBInfo(UniProt.ws::UniProt.ws(), extIDs, col = "xref_kegg")
  names(genesfupkb) = c("extid", "keggid")
  upids = allextIDs[which(allextIDs$extid %in% genesfupkb$extid),]
  genesfupkb =  merge.data.frame(genesfupkb, upids, by = "extid")
  names(genesfupkb) = c("uniprotid", "keggid", "internalid", "source")
  return(genesfupkb)
}

