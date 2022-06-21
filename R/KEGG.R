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
  if ("NCBI-GeneID" %in% colnames(genes))
    genes = genes[,-which(colnames(genes) %in% c("NCBI-GeneID"))]
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
#' @param geneList A list of KEGG genes to add
#'
#' @return A Mully graph with the all layers,nodes and edges
#' @export
#'
#' @examples
#' g= mully("Test")   #This is a new graph but it is allowed to have an existing graph with nodes and edges already added
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' up = UniProt.ws()
#' graph=addGenesLayer(g,geneList)
addGenesLayer <- function(g, geneList) {
  g = addLayer(g, "KEGGGenes")
  genes = getKeggGene(geneList)
  genes = interactionKegg(genes)
  genes = simplifyInteractionKegg(genes)
  
  nbnodes = 0 #count of nodes added
  nbedges = 0 #count of edges added
  #Add Genes' Nodes
  message("Multipath: Adding GENE Nodes")
  nodesToAdd=unique(genes[,1])
  for(i in 1:length(nodesToAdd)){
    if (is.null (getNode(g, nodesToAdd[i]))) {
      nbnodes = nbnodes + 1
      g = mully::addNode(
        g,
        nodeName = paste0("KEGGGene",nbnodes),
        layerName = "KEGGGenes",
        # attributes = genes[i,4],
        attributes =  list("id"=nodesToAdd[i],"database" = "KEGG")
      )
          message("Multipath: DONE - GENE Nodes Added")
    }
  }   
  
  message("Added Nodes = ", nbnodes)
  return(g)
}
# addAllGenesEdges <- function(g,geneList){
# 
# internalId= getInternalIDs(g)
# message("Multipath: Adding Edges")
# for (i in 1:length(genes[,1])) {
#   startName= internalId[which(internalId$extid==genes[i,1]),]$id
#   endname = internalId[which(internalId$extid==genes[i,2]),]$id
#   for (j in 1:length(endname)){
#     endName=endname[j]
#     print(paste0(startName,endName,sep=" and "))
#     if(!is.na(endName) & !is.null(startName) & !is.null(endName) & length(startName)!=0 & length(endName)!=0 & getIDEdge(g,startName,endName) == 0 ){
#       g = mully::addEdge(g,nodeStart=startName,nodeDest=endName,list("source" = "kegggenes"))
#       nbedges = nbedges + 1
#       message("Multipath: DONE - GENE Edge Added")
#       
#     }
#   }
# }

#  
#   message(" Added Edges = ", nbedges)
#   return(g)
# }
addGenesEdges <- function(g,relation){
  nbedges = 0 #count of edges added
  internalId= getInternalIDs(g)
  message("Multipath: Adding Edges")
  for (i in 1:length(relation[,1])) {
    startName= V(g)[which(V(g)$id==relation[i,1])]$name
    endname = internalId[which(internalId$extid==relation[i,2]),2]
    for (j in 1:length(endname)){
      endName=endname[j]
      if(!is.na(endName) && !is.null(startName) & !is.null(endName) & length(startName)!=0 & length(endName)!=0 ){
        g = mully::addEdge(g,nodeStart=startName,nodeDest=endName,list("source" = relation[i,3] ))
        nbedges = nbedges + 1
        message("Multipath: DONE - GENE Edge Added")
        
      }
    }
  }
  
  message(" Added Edges = ", nbedges)
  return(g)
}


#' Get KEGG Genes to UniProt relations from KEGG
#'
#' @param geneList    List of Genes as formatted in KEGG Genes 
#' @param proteinList List of Proteins as formatted in UPKB
#'
#' @return Dataframe of the Genes-Proteins relation
#' @export
#' @importFrom dplyr select
#'
#' @examples
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' kEGG2UPKB = getKEGGtoUPKB(geneList,proteinList)

getKEGGtoUPKB <- function(geneList, proteinList) {
  genesKegg = getKeggGene(geneList)
  genesKegg = interactionKegg(genesKegg)
  genesKegg = simplifyInteractionKegg(genesKegg)
  genes = (dplyr::filter(genesKegg, dbName == "UniProt"))
  allRelations = dplyr::select(genes, c("keggEntry", "dbId","source"))
  names(allRelations) = c("keggid", "upid","source")
  allRelations = allRelations[order(allRelations$upid), ]
  if (missing(geneList)) {
    return(allRelations)
  }
  relations = allRelations[which(allRelations$'upid' %in% proteinList), ]
  return(relations)
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
getUPKBtoKEGG<-function(up,proteinList,geneList){
  allRelations=getUPKBInfo(up,proteinList,col = c("UNIPROTKB","KEGG"))
  allRelations = cbind(allRelations,"source"="UPKB")
  names(allRelations)=c("upid","keggid","source")
  allRelations=allRelations[order(allRelations$upid),]
  if(missing(geneList)){
    return(allRelations)
  }
  relations=allRelations[which(allRelations$'keggid'%in%geneList),]
  return(relations)
}

#' Get KEGG Genes to UniProt relations from KEGG
#'
#' @param geneList    List of Genes as formatted in KEGG Genes 
#' @param proteinList List of Proteins as formatted in UPKB
#'
#' @return Dataframe of the Genes-Proteins relation
#' @export
#'
#' @examples
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' kEGG2UPKB = getKEGGtoUPKB(geneList,proteinList)

getKEGGtOMIM<- function(geneList, proteinList) {
  genesKegg = getKeggGene(geneList)
  genesKegg = interactionKegg(genesKegg)
  genesKegg = simplifyInteractionKegg(genesKegg)
  genes = (filter(genesKegg, dbName == "UniProt"))
  allRelations = dplyr::select(genes, c("keggEntry", "dbId"))
  names(allRelations) = c("keggid", "upid")
  allRelations = allRelations[order(allRelations$upid), ]
  if (missing(geneList)) {
    return(allRelations)
  }
  relations = allRelations[which(allRelations$'upid' %in% proteinList), ]
  return(relations)
}
# library(biodb)
# library(classyfireR)

#get other database cross reference when found from upkb by dbABBRV
getUPKBtoCrossReference <- function(proteinList, dbAbbrv) {
  df_upkbtoany = data.frame(
    Entry = is.character(c()),
    CrossReference = is.character(c()),
    stringsAsFactors = FALSE
  )[-1, ]
  h = 1
  i = 1
  while (i <= length(proteinList)) {
    #Set of 100 elements
    indexEnd = h * 100
    if (length(proteinList) <= indexEnd)
      indexEnd = length(proteinList)
    #Check if out of bound
    if (indexEnd > length(proteinList))
      indexEnd = length(proteinList)
    list_protein =  c(proteinList[i:indexEnd])
    df_queryentry = queryfromUPKB(list_protein, dbAbbrv)
    df_upkbtoany = merge.data.frame(
      x = df_upkbtoany,
      y = df_queryentry,
      all.x = TRUE,
      all.y = TRUE
    )
    i = indexEnd + 1
    h = h + 1
  }
  df_upkbtoany = df_upkbtoany[,-which(colnames(df_upkbtoany) %in% c("CrossReference"))]
  return(df_upkbtoany)
}
queryfromUPKB <- function(proteinList, dbAbbrv) {
  queryentry = c()
  queryentrybtwn = "+or+"
  dbname = dbAbbrv
  for (i in 1:length(proteinList)) {
    queryentrytemp = noquote(proteinList[i])
    queryentrytemp = paste(noquote("id:"), queryentrytemp, sep = "")
    if (i == length(proteinList))
      queryentry = paste(queryentry, queryentrytemp, sep = "")
    if (i < length(proteinList))
      queryentry = paste(queryentry, queryentrytemp, queryentrybtwn, sep = "")
  }
  generatelink =  paste (
    noquote("https://www.uniprot.org/uniprot/?query="),
    queryentry,
    noquote("&format=tab&columns=id,"),
    noquote("database("),
    dbname,
    noquote(")"),
    sep = ""
  )
  result = fread(generatelink)
  df_queryentry = result
  return(df_queryentry)
}

getKEGGGenesUPKBRelations <- function(up, proteinList, geneList) {
  upkbtokegg = na.omit(getUPKBtoKEGG(up, proteinList, geneList))
  upkbtokegg = upkbtokegg[-which(upkbtokegg[,1] == "NA" | upkbtokegg[,2] == "NA"),]
  keggtoupkb = na.omit(getKEGGtoUPKB(geneList, proteinList))
  keggtoupkb = keggtoupkb[-which(keggtoupkb[,1] == "NA" | keggtoupkb[,2] == "NA"),]
  relations = merge.data.frame(upkbtokegg,
                               keggtoupkb,
                               by = c("keggid", "upid"),
                               all = T)
  for (i in 1:length(relations$keggid)) {
    relations$source[i] = ""
    if (!is.na(relations$'source.x'[i]))
      relations$source[i] = paste(relations$'source.x'[i], " ")
    if (!is.na(relations$'source.y'[i]))
      relations$source[i] = paste(relations$source[i], relations$'source.y'[i])
    relations$source[i] = str_trim(relations$source[i], side = "right")
  }
  relations = relations[, c("keggid", "upid", "source")]
  relations = relations[order(relations$keggid),]
  return(relations)
}

getUPKBDBRelations<-function(up,data,proteinList,drugList){
  upkbtodb=getUPKBtoDB(up,proteinList,drugList)
  dbtoupkb=getDBtoUPKB(data,drugList,proteinList)
  lnames=intersect(names(upkbtodb),names(dbtoupkb))
  upkbtodb$source="UniProt"
  dbtoupkb$source="DrugBank"
  relations=merge.data.frame(upkbtodb,dbtoupkb[ ,c("dbid","upid","type","dbproteinid","source")],by=lnames,all=T)
  relations$source=paste(relations$'source.x'," ",relations$'source.y')
  relations=relations[,c("dbid","upid","type","dbproteinid","source")]
  relations=relations[order(relations$dbid),]
  return(relations)
}

getRelatedGenes <- function (g,biopax,layer,database){
  allextIDs = getExternalIDs(biopax,getLayer(g,layer)$name)
  extIDs = unique(allextIDs[which(allextIDs$database==database),]$extid)
  genesfupkb=getUPKBInfo(UniProt.ws(),extIDs,col = c("UNIPROTKB","KEGG"))
  
  return(genesfupkb)
}