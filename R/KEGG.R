# Get A Gene using KeggGet
#' Get Kegg gene
#'
#' @param geneList  List of Genes as formatted in KEGG Genes
#'
#' @return A data frame of kegg data
#' @export
#'
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' genes=  getKeggGene(geneList)
#'
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
#'
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
#'
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
#'
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
#'
#' @examples
#' g= mully("Test")   #This is a new graph but it is allowed to have an existing graph with nodes and edges already added
#' biopax= readBiopax("TCR.R-HSA-202403_level2.owl") #The file is parsed using rbiopaxparser
#' graph=addGenesLayer(g,geneList)

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
  print(startName)
  endName = data$kegginternalid
  print(endName)
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

#' @param geneList    List of Genes as formatted in KEGG Gene

#' Get KEGG Genes to UniProt relations from KEGG
#'s
#' @param proteinList List of Proteins as formatted in UPKB
#' @param dbName      Name of database to find as string- only 1 input allowed
#'
#' @return Dataframe of the Genes-Proteins relation
#' @export
#' @importFrom dplyr select
#'
#' @examples
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' kEGG2UPKB = getKEGGtoDATABASE(geneList,proteinList)

getKEGGtoDATABASE <-
  function(dbName, geneList, toDatabaseList = optional) {
    genesKegg = getKeggGene(geneList)
    genesKegg = interactionKegg(genesKegg)
    genesKegg = simplifyInteractionKegg(genesKegg)
    genes = applyFilter(genesKegg, dbName)
    allRelations = dplyr::select(genes, c("keggEntry", "dbId", "source", "attributes"))
    names(allRelations) = c("keggid", dbName, "source", "keggattributes")
    allRelations = allRelations[order(allRelations[[dbName]]), ]
    if (missing(toDatabaseList)) {
      return(allRelations)
    }
    relations = allRelations[which(allRelations[[dbName]] %in% toDatabaseList), ]
    return(relations)
  }

applyFilter <- function(dataFrame, columnName) {
  if (tolower(columnName) == "uniprot")
    df = (dplyr::filter(dataFrame, dbName == "UniProt"))
  if (tolower(columnName) == "omim")
    df = (dplyr::filter(dataFrame, dbName == "OMIM"))
  if (tolower(columnName) == "ensembel")
    df = (dplyr::filter(dataFrame, dbName == "Ensembl"))
  
  return(df)
}

#' Get UniProt Proteins to KEGG Genes relations from UPKB
#'
#' @param up          The UniProt.ws Object
#' @param proteinList List of Proteins as formatted in UPKB
#' @param geneList    List of Genes as formatted in KEGG Genes
#'
#' @return Dataframe of the Proteins-Genes relation
#' @export
#'
#' @examples
#' up = UniProt.ws()
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' UPKB2KEGG = getUPKBtoKEGG(up,geneList,proteinList)
#'
getUPKBtoKEGG <- function(up, proteinList, geneList) {
  allRelations = getUPKBInfo(up, proteinList, col = "xref_kegg")
  allRelations = cbind(allRelations, "source" = "UPKB")
  names(allRelations) = c("UniProt", "keggid", "source")
  allRelations = allRelations[order(allRelations$UniProt),]
  if (missing(geneList)) {
    return(allRelations)
  }
  relations = allRelations[which(allRelations$'keggid' %in% geneList),]
  return(relations)
}

#' Gets the genes and proteins that are referenced to each other
#'
#' @param g            A Mully graph object
#' @param biopax   A biopax object
#'
#' @return              Dataframe of related proteins and genes
#' @export
#'
#' @examples
getKeggUpkbRelations <- function(g, biopax) {
  upkbtokeggDf = getRelatedGenes(g, biopax)
  geneList = upkbtokeggDf$keggid
  keggtoupkbDf = getKEGGtoDATABASE("UniProt", geneList)
  names(keggtoupkbDf) = c("keggid", "uniprotid", "source", "keggattributes")
  relations = merge.data.frame(upkbtokeggDf,
                               keggtoupkbDf,
                               by = c("keggid", "uniprotid"),
                               all = T)
  relations = na.omit(relations)
  return(relations)
}

#' Gets the genes and proteins that are referenced to each other
#'
#' @param up            The UniProt.ws Object
#' @param proteinList   List of Proteins as formatted in UPKB
#' @param geneList      List of Genes as formatted in KEGG Genes
#'
#' @return              Dataframe of related proteins and genes
#' @export
#'
#' @examples
getKEGGGenesUPKBRelations <- function(up, proteinList, geneList) {
  upkbtokegg = getUPKBtoKEGG(up, proteinList, geneList)
  #upkbtokegg = upkbtokegg[-which(upkbtokegg[, 1] == "NA" |upkbtokegg[, 2] == "NA"), ]
  
  keggtoupkb = getKEGGtoDATABASE("UniProt", geneList, proteinList)
  #keggtoupkb = keggtoupkb[-which(keggtoupkb[, 1] == "NA" | keggtoupkb[, 2] == "NA"), ]
  
  relations = merge.data.frame(upkbtokegg,
                               keggtoupkb,
                               by = c("keggid", "UniProt"),
                               all = T)
  for (i in 1:length(relations$keggid)) {
    relations$source[i] = ""
    if (!is.na(relations$'source.x'[i]))
      relations$source[i] = paste(relations$'source.x'[i], " ")
    if (!is.na(relations$'source.y'[i]))
      relations$source[i] = paste(relations$source[i], relations$'source.y'[i])
    relations$source[i] = str_trim(relations$source[i], side = "right")
  }
  relations = relations[, c("keggid", "UniProt", "source")]
  relations = relations[order(relations$keggid),]
  return(relations)
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

#' Get KEGG Genes to Omim relations from KEGG
#'
#' @param geneList    List of Genes as formatted in KEGG Genes
#'
#' @return Dataframe of the Genes-Disease relation
#' @export
#'
#' @examples
#'
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' kEGG2OMIM = getKEGGtoOMIM(geneList,proteinList)

getKEGGtoOMIM <- function(geneList) {
  genesKegg = getKeggGene(geneList)
  genesKegg = interactionKegg(genesKegg)
  genesKegg = simplifyInteractionKegg(genesKegg)
  genes = dplyr::filter(genesKegg, dbName == "OMIM")
  genesToOmim = dplyr::select(genes, c("keggEntry", "dbId", "source"))
  return(genesToOmim)
}

#' Get Omim to KEGG Genes relations from OMIM
#'
#' @param omimIds list of OMIM Ids
#'
#' @return Dataframe of the Genes-OMIM relation
#' @export
#'
#'@Note should be preceded by calling  set_key('KEY')
#' @examples
getOmimToKEGG <- function(omimIds) {
  df_omim = data.frame(
    OMIM = c(),
    KEGGID = c(),
    source = c(),
    attributes = c(),
    stringsAsFactors = FALSE
  )
  for (i in 1:length(omimIds)) {
    entry = as.character(omimIds[i])
    omim_result = get_omim(entry, geneMap = TRUE)
    omim_data = xmlToList(omim_result)
    omim_title = omim_data$entryList$entry$titles$preferredTitle
    omim_gene = omim_data$entryList$entry$geneMap$geneIDs
    translateGeneID2KEGGID = translateGeneID2KEGGID(omim_gene)
    new_row = data.frame(
      OMIM = entry,
      KEGGID = translateGeneID2KEGGID,
      source = "OMIM",
      attributes = omim_title
    )
    df_omim = rbind(df_omim, new_row)
  }
  return(df_omim)
}

#' Gets the genes and diseases that are referenced to each other. Genes are extracted from the graph
#'
#' @param g A Mully Graph Object
#' @param biopax A Biopax Object
#'
#' @return A Dataframe with kegggenes ids and omim ids
#' @export
#'
#' @note must be preceded by addGenesLayer(g,biopax) function
#'
#' @examples
#' biopax=readBiopax("wnt.owl")
#' pathwayID=listPathways(biopax)$id[1]
#' g=Multipath::pathway2Mully(biopax,pathwayID)
#' g=addGenesLayer(g,biopax)
#' dataframe = getKeggOmimRelation(g,biopax)
getKeggOmimRelation <- function(g, biopax) {
  gene_list = getLayer(g, "KEGGGenes")$id
  gene_internalid = getLayer(g, "KEGGGenes")$name
  df_gene_internalid = data.frame(KeggEntry = gene_list,
                                  gene_internalid = gene_internalid)
  keggToOmim = getKEGGtoOMIM(gene_list)
  df_gene_internalid = df_gene_internalid[which(keggToOmim$keggEntry  %in%  my_data$KeggEntry), ]$gene_internalid
  keggToOmim = cbind(keggToOmim, df_gene_internalid)
  names(keggToOmim) = c("KEGGID", "OMIM", "source", "geneInternalId")
  getOmimToKEGG = getOmimToKEGG(keggToOmim$OMIM)
  relations = merge.data.frame(keggToOmim,
                               getOmimToKEGG,
                               by = c("KEGGID", "OMIM"),
                               all = T)
  relations$source = paste(relations$source.x, relations$source.y, sep = "; ")
  relations = relations %>% select(-source.x, -source.y)
  
  relations = na.omit(relations)
  return(relations)
}

#' Get Omim to UPKB relations from OMIM
#'
#' @param omimIds list of OMIM Ids
#'
#' @return Dataframe of the Proteins-OMIM relation
#' @export
#'
#'@Note should be preceded by calling  set_key('KEY')
#' @examples
getOmimToUPKB <- function(omimIds) {
  df_omim = data.frame(
    OMIM = c(),
    UPid = c(),
    source = c(),
    attributes = c(),
    stringsAsFactors = FALSE
  )
  for (i in 1:length(omimIds)) {
    entry = as.character(omimIds[i])
    omim_result = get_omim(entry, geneMap = TRUE)
    omim_data = xmlToList(omim_result)
    omim_title = omim_data$entryList$entry$titles$preferredTitle
    omim_protein  = omimdata$entryList$entry$externalLinks$swissProtIDs
    new_row = data.frame(
      OMIM = entry,
      Upid = omim_protein,
      source = "OMIM",
      attributes = omim_title
    )
    df_omim = rbind(df_omim, new_row)
  }
  return(df_omim)
}


#' Add OMIM Layer with its respective Nodes
#'
#' @param g A Mully graph
#' @param biopax A Biopax object
#'
#' @return g A Mully graph including OMIM nodes and possible edges
#' @export
#'
#' @examples
#' g= mully("Test")   #This is a new graph but it is allowed to have an existing graph with nodes and edges already added
#' biopax= readBiopax("TCR.R-HSA-202403_level2.owl") #The file is parsed using rbiopaxparser
#' g=addGenesLayer(g,biopax)
#' g=addDiseaseLayer(g,biopax)
addDiseaseLayer <- function(g, biopax) {
  g = addLayer(g, "OMIM")
  if (!isLayer(g, "KEGGGenes"))
    g = addGenesLayer(g, biopax)
  
  data = getKeggOmimRelation(g, biopax)
  omimIds = data$OMIM
  keggIds = data$KEGGID
  source = data$source
  attributes = data$attributes
  nbnodes = 0 #count of nodes added
  nbedges = 0 #count of edges added
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
          "database" = source[i],
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
  #Adding OMIM-Proteins Edges
  
  nbedges = 0 #count of edges added
  df_omimfromUPKB = getUPKBRelatedDiseases(g, biopax)
  for (i in 1:length(df_omimfromUPKB$omimid)) {
    nodeName = which(getLayer(g, "OMIM")$id == df_omimfromUPKB$omimid[i])
    print(nodeName)
    if (length(nodeName) == 0) {
      nbnodes = nbnodes + 1
      g = addNode(
        g,
        nodeName = paste0("OMIM", nbnodes),
        layerName = "OMIM",
        attributes = list("id" = df_omimfromUPKB$omimid[i])
      )
      print("ADDED")
    }
  }
  
  for (i in 1:nrow(df_omimfromUPKB)) {
    omim_internalid = getLayer(g, "OMIM")$name[which(getLayer(g, "OMIM")$id == df_omimfromUPKB$omimid[i])]
    protein_internalid = df_omimfromUPKB$protein_internalid[i]
    g = mully::addEdge(g,
                       nodeStart = omim_internalid,
                       nodeDest = protein_internalid,
                       list("source" = "UniProt;OMIM"))
    nbedges = nbedges + 1
    message(
      "Multipath: DONE - Disease Edge Added  ",
      paste(omim_internalid, protein_internalid, sep = "<>")
    )
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
#'
#' @examples
getUPKBRelatedDiseases <- function (g, biopax) {
  allextIDs = getExternalIDs(biopax, getLayer(g, "protein")$name)
  extIDs = unique(allextIDs[which(allextIDs$database == "UniProt"),]$extid)
  omimfupkb = getUPKBInfo(UniProt.ws::UniProt.ws(), extIDs, col = "xref_mim")
  names(omimfupkb) = c("extid", "omimid")
  df_omimfupkb = data.frame(extid = c(),
                            omimid = c())
  
  for (i in 1:nrow(omimfupkb)) {
    omimid = omimfupkb$omimid[i]
    
    # split the omimid into split of six characters
    split = substring(omimid, seq(1, nchar(omimid), 6),
                      seq(6, nchar(omimid) + 5, 6))
    for (i in 1:length(split)) {
      new_row = data.frame(extid = omimfupkb$extid[i], omimid = split[i])
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
