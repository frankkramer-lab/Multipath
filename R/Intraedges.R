#' Get UniProt Proteins to DrugBank Drugs relations from UniProt
#'
#' @param proteinList The list of UniProt Ids of the proteins
#' @param drugList The ID of the DrugBank drug entry starting with "DB". This argument can be either a string (one drug) or a list of strings (multiple drugs).
#'
#' @return A dataframe containing the connections between UniProt proteins and DrugBank drugs retrieved from UniProt
#' @export
#'
#' 
#' @examples
#' \dontrun{ 
#' getUPKBtoDB(c("P02747","P00734","P07204"),c("DB00001","DB00002"))
#' }
getUPKBtoDB<-function(proteinList,drugList){
  allRelations=getUPKBInfo(proteinList,col = c("xref_drugbank"))
  result=data.frame(upid=is.character(c()),dbid=is.character(c()),stringsAsFactors = F)[-1,]
  names(allRelations)=c("upid","dbid")
  allRelations=allRelations[order(allRelations$upid),]
  rows=dim(allRelations)[1]
  for(i in 1:rows){
    listInter=as.list(str_split(allRelations[i,2],";"))[[1]]
    for(j in 1:length(listInter)){
      if("" == listInter[j])
        next()
      entry=c(allRelations[i,1],listInter[j])
      result[dim(result)[1]+1,]=entry
    }
  }
  if(missing(drugList)){
    return(result)
  }
  result=result[which(result$'dbid'%in%drugList),]
  return(result)
}

#' Get DrugBank Drugs to UniProt Proteins Relations from DrugBank
#'
#' @param data  The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#' @param proteinList The list of UniProt Ids of the proteins
#'
#' @return A dataframe containing the connections between DrugBank drugs and UniProt proteins retrieved from DrugBank
#' @export
#'
#' @examples
#' \dontrun{ 
#' data=readDBXML(DBXMLFilePath)
#' getDBtoUPKB(data,c("DB00001","DB00002","DB00006"),c("P02747","P00734","P07204","P05164"))
#' }
getDBtoUPKB<-function(data,drugList,proteinList){
  l=list()
  
  #Get Targets
  allTargets=getDBTargets(data,drugList)
  targets=allTargets[which(allTargets$'upid'%in%proteinList),]
  if(dim(targets)[1]>0){
    targets$type="Target"
    l[[length(l)+1]]=targets
  }
  
  #Get Enzymes
  allEnzymes=getDBEnzymes(data,drugList)
  enzymes=allEnzymes[which(allEnzymes$'upid'%in%proteinList),]
  if(dim(enzymes)[1]>0){
    enzymes$type="Enzyme"
    l[[length(l)+1]]=enzymes
  }
  
  #Get Transporters
  allTransporters=getDBTransporters(data,drugList)
  transporters=allTransporters[which(allTransporters$'upid'%in%proteinList),]
  if(dim(transporters)[1]>0){
    transporters$type="Transporter"
    l[[length(l)+1]]=transporters
  }
  
  #Get Carriers
  allCarriers=getDBCarriers(data,drugList)
  carriers=allCarriers[which(allCarriers$'upid'%in%proteinList),]
  if(dim(carriers)[1]>0){
    carriers$type="Carrier"
    l[[length(l)+1]]=carriers
  }
  
  n=length(l)
  if(n==0){
    warning("No relations between the given nodes")
    return(NULL)
  }
  relations=l[[1]]
  if(n==1)
    return(l[[1]])
  for(i in 2:n){
   lnames=intersect(names(relations),names(l[[i]]))
   relations=merge.data.frame(relations,l[[i]],by=lnames,all=T) 
  }
  relations=relations[order(relations$dbid),]
  return(relations)
}

#' Get Protein and Drugs relations from UniProt and DrugBank
#'
#' @param data The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param proteinList The list of UniProt Ids of the proteins
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#'
#' @return A dataframe containing the connections between DrugBank drugs and UniProt proteins retrieved from DrugBank and UniProt
#' @export
#' 
#' 
#' @note
#' Should be preceded by:
#' loadDBXML(DrugBankFile) to get the argument data
#' @examples
#' \dontrun{
#' data=readDBXML(DBXMLFilePath)
#' relations=getUPKBDBRelations(data,c("P02747","P05164"),c("DB00001","DB00006"))
#' }
getUPKBDBRelations<-function(data,proteinList,drugList){
  upkbtodb=getUPKBtoDB(proteinList,drugList)
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

#' Get UniProt Proteins to KEGG Genes relations from UPKB
#'
#' @param proteinList List of Proteins as formatted in UPKB
#' @param geneList    List of Genes as formatted in KEGG Genes
#'
#' @return Dataframe of the Proteins-Genes relation
#' @export
#' @author Mohammad Al Maaz, Zaynab Hammoud
#' @examples
#' proteinList = c("P02747","P00734","P07204","A0A0S2Z4R0","O15169")
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' UPKB2KEGG = getUPKBtoKEGG(geneList,proteinList)
getUPKBtoKEGG <- function(proteinList, geneList) {
  allRelations = getUPKBInfo(proteinList, col = "xref_kegg")
  allRelations = cbind(allRelations, "source" = "UPKB")
  names(allRelations) = c("UniProt", "keggid", "source")
  allRelations = allRelations[order(allRelations$UniProt),]
  if (missing(geneList)) {
    return(allRelations)
  }
  allRelations$'keggid'=gsub(";", "",allRelations$'keggid')
  relations = allRelations[which(allRelations$'keggid' %in% geneList),]
  return(relations)
}

#' Get KEGG Genes to a Database cross-reference  from KEGG
#'
#' @param geneList List of Proteins as formatted in UPKB
#' @param dbName      Name of database to find as string- only 1 input allowed
#'
#' @return Dataframe of the Genes-Proteins relation
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' keggToDb = getKEGGtoDATABASE("UniProt",geneList)
getKEGGtoDATABASE <-
  function(dbName, geneList) {
    genesKegg = getKeggGene(geneList)
    genesKegg = interactionKegg(genesKegg)
    genesKegg = simplifyInteractionKegg(genesKegg)
    if (tolower(dbName) == "uniprot")
      genes = (dplyr::filter(genesKegg, dbName == "UniProt"))
    if (tolower(dbName) == "omim")
      genes = (dplyr::filter(genesKegg, dbName == "OMIM"))
    if (tolower(dbName) == "ensembel")
      genes = (dplyr::filter(genesKegg, dbName == "Ensembl"))
    allRelations = dplyr::select(genes, c("keggEntry", "dbId", "source", "attributes"))
    names(allRelations) = c("keggid", dbName, "source", "keggattributes")
    allRelations = allRelations[order(allRelations[[dbName]]), ]
    return(allRelations)
  }

#' Get KEGG Genes to Omim relations from KEGG
#'
#' @param geneList    List of Genes as formatted in KEGG Genes
#'
#' @return Dataframe of the Genes-Disease relation
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' geneList=c("hsa:122706","hsa:4221","hsa:8312")
#' kEGG2OMIM = getKEGGtoOMIM(geneList)
getKEGGtoOMIM <- function(geneList) {
  genesKegg = getKeggGene(geneList)
  genesKegg = interactionKegg(genesKegg)
  genesKegg = simplifyInteractionKegg(genesKegg)
  genes1 = genesKegg[which(genesKegg$dbName=="OMIM"),]
  genesToOmim =
    dplyr::select(genes1, c("keggEntry", "dbId", "source"))
  return(genesToOmim)
  
}

#' Gets the genes and proteins that are referenced to each other
#'
#' @param g            A Mully graph object
#' @param biopax   A biopax object
#'
#' @return              Dataframe of related proteins and genes
#' @export
#' @author Mohammad Al Maaz
#' @examples
#' \dontrun{ 
#' biopax=readBiopax("wnt.owl") 
#' pathwayID=listPathways(biopax)$id[1]
#' g=Multipath::pathway2Mully(biopax,pathwayID) #This is a mully graph that contains a protein layer
#' g=addGenesLayer(g,biopax)
#' UPKB2KEGG = getUPKBtoKEGG(g, biopax)
#' }
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

#' Get Omim to UPKB relations from OMIM
#'
#' @param omimIds list of OMIM Ids
#'
#' @return Dataframe of the Proteins-OMIM relation
#' @export
#'
#' @note  should be preceded by calling  romim::set_key('KEY'). The KEY could be requested via omim's official website.
#' @examples
#' \dontrun{   
#' library(romim)
#' romim::set_key('KEY')
#' OmimToUPKB = getOmimToUPKB(c("611137", "613733", "603816"))
#' }
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
    omim_result = get_omim(entry, externalLinks = TRUE)
    omim_data = XML::xmlToList(omim_result)
    omim_title = omim_data$entryList$entry$titles$preferredTitle
    omim_protein  = omim_data$entryList$entry$externalLinks$swissProtIDs
    
    if (!is.null(omim_protein)) {
      omim_protein = unlist(str_split(omim_protein, ","))
      for (i in 1:length(omim_protein)) {
        new_row = data.frame(
          OMIM = entry,
          UPid = omim_protein[i],
          source = "OMIM",
          attributes = omim_title
        )
        
        df_omim = rbind(df_omim, new_row)
      }
    }
  }
  return(df_omim)
}

#' Get Omim to KEGG Genes relations from OMIM
#'
#' @param omimIds list of OMIM Ids
#' @importFrom romim get_omim
#' @return Dataframe of the Genes-OMIM relation
#' @export
#' @importFrom XML xmlToList
#' @importFrom  xml2 read_xml
#' @importFrom  XML xmlParse
#'
#' @note  should be preceded by calling  romim::set_key('KEY'). The KEY could be requested via omim's official website.
#' @author Mohammad Al Maaz
#'
#' @examples 
#' \dontrun{ 
#' library(romim)
#' romim::set_key('KEY')
#' kEGG2OMIM = getOmimToKEGG(c("611137", "613733", "603816"))
#' }
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
    omim_result = romim::get_omim(entry, geneMap = TRUE)
    omim_data = xmlToList(omim_result)
    omim_title = omim_data$entryList$entry$titles$preferredTitle
    omim_gene = omim_data$entryList$entry$geneMap$geneIDs
    translateGeneID2KEGGID = KEGGgraph::translateGeneID2KEGGID(omim_gene)
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
#' @author Mohammad Al Maaz
#' @examples
#' \dontrun{ 
#' biopax=readBiopax("wnt.owl")
#' pathwayID=listPathways(biopax)$id[1]
#' g=Multipath::pathway2Mully(biopax,pathwayID)
#' g=addGenesLayer(g,biopax)
#' dataframe = getKeggOmimRelation(g,biopax)
#' }
getKeggOmimRelation <- function(g, biopax) {
  gene_list = getLayer(g, "KEGGGenes")$id
  gene_internalid = getLayer(g, "KEGGGenes")$name
  df_gene_internalid = data.frame(KeggEntry = gene_list,
                                  gene_internalid = gene_internalid)
  keggToOmim = getKEGGtoOMIM(gene_list)
  df_gene_internalid = df_gene_internalid[which(keggToOmim$keggEntry  %in%  df_gene_internalid$KeggEntry), ]$gene_internalid
  keggToOmim = cbind(keggToOmim, df_gene_internalid)
  names(keggToOmim) = c("KEGGID", "OMIM", "source", "geneInternalId")
  View(keggToOmim)
  getOmimToKEGG = getOmimToKEGG(keggToOmim$OMIM)
  relations = merge.data.frame(keggToOmim,
                               getOmimToKEGG,
                               by = c("KEGGID", "OMIM"),
                               all = T)
  relations$source = paste(relations$source.x, relations$source.y, sep = "; ")
  relations = relations %>% dplyr::select(-source.x, -source.y)
  
  relations = na.omit(relations)
  return(relations)
}
