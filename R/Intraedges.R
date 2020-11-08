#' Get UniProt Proteins to DrugBank Drugs relations from UniProt
#'
#' @param up The UniProt.ws Object
#' @param proteinList The list of UniProt Ids of the proteins
#' @param drugList The ID of the DrugBank drug entry starting with "DB". This argument can be either a string (one drug) or a list of strings (multiple drugs).
#'
#' @return A dataframe containing the connections between UniProt proteins and DrugBank drugs retrieved from UniProt
#' @export
#'
#' 
#' @note
#' Should be preceded by UniProt.ws() to get the UniProt.ws Object
#' @examples
#' \dontrun{ 
#' up=UniProt.ws()
#' getUPKBtoDB(up,c("P02747","P00734","P07204"),c("DB00001","DB00002"))
#' }
getUPKBtoDB<-function(up,proteinList,drugList){
  allRelations=getUPKBInfo(up,proteinList,col = c("UNIPROTKB","DRUGBANK"))
  relations=allRelations[which(allRelations$'DRUGBANK'%in%drugList),]
  names(relations)=c("upid","dbid")
  relations=relations[order(relations$upid),]
  return(relations)
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
#' @param up The UniProt.ws Object
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
#' 1. UniProt.ws() to get the UniProt.ws Object
#' 2. loadDBXML(DrugBankFile) to get the argument data
#' @examples
#' \dontrun{
#' up=UniProt.ws() 
#' data=readDBXML(DBXMLFilePath)
#' relations=getUPKBDBRelations(up,data,c("P02747","P05164"),c("DB00001","DB00006"))
#' }
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