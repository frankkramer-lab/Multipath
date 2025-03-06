#' Load DrugBank XML file
#' 
#' @param file The path to the DrugBank XML file. This can be downloaded from the DrugBank official Website (drugbank.ca). An account with an institutional e-mail is required.
#'
#' @return A dataframe containing the parsed information from DrugBank. This can be used to extract any additional information on the DrugBank entries
#' @export
#'
#' @note
#' This function should be called before using any function to query the DrugBank database. Since the parsing of DrugBank takes time, this function should only be called once.
loadDBXML<-function(file){
  vectoroptions=dbparser::drug_node_options()
  data=dbparser:::parseDrugBank(file, drug_options =vectoroptions)
  return(data$drugs)
}

#' Get DrugBank drug entry
#'
#' @param data The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile) 
#' @param drug The ID or list of IDs of the DrugBank drug entries starting with "DB"
#'
#' @return A dataframe containing the DrugBank entry with its information
#' @export
#'
#' @examples
#' \dontrun{ 
#' data=loadDBXML(DrugBankFilePath)
#' getDBDrug(data, "DB00001")
#' }
getDBDrug<-function(data,drug){
  drugs=as.data.frame(data[[1]])
  entry=drugs[which(drugs$'drugbank_id'%in%drug),]
  if(length(entry)==0){
    warning("The given DrugBank ID does not exist")
    return(NULL)
  }
  if(dim(entry)[1]!=length(drug))
    warning("Some of the given IDs are not DrugBank IDs and will be skipped.")
  names(entry)[which(names(entry)=="name")]="dbname"
  return(entry)
}

#' Get DrugBank Drug to Drug Interactions 
#' 
#' @param data The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile) 
#' @param drug The ID of the DrugBank drug entry starting with "DB". This argument can be either a string (one drug) or a list of strings (multiple drugs).
#'
#' @return A dataframe containing the DrugBank interactions in which the given drug is involved
#' @export
#'
#' @examples 
#' \dontrun{ 
#' data=loadDBXML(DBXMLFilePath)
#' getDBDrugInteractions(data,"DB06605")
#' }
getDBDrugInteractions<-function(data,drug){
  if(!is.list(drug))
    drugs=c(drug)
  #remove duplicate entries
  drugs=unique(drugs)
  #remove non-drug entries
  checkEntries=getDBDrug(data,drugs)
  isDBEntry=drugs%in%checkEntries[,1]
  nonDrugs=drugs[isDBEntry==F]
  drugs=drugs[isDBEntry==T]
  if(length(drugs)==0)
    stop("The given drug IDs are false.")
  if(length(nonDrugs)!=0)
    warning(paste("The following arguments are not DrugBank IDs and will be skipped:"),nonDrugs)
  #Get the interactions
  interactions=as.data.frame(data$'drug_interactions')
  drugInteractions=interactions[which(interactions$'drugbank_id'%in%drugs | interactions$'target_drugbank_id'%in%drugs),]
  return(drugInteractions)
}

#' Add a drug layer to a mully graph
#'
#' @param g The mully graph
#' @param drugList The list of DrugBank Ids of the drugs to be added. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#' @param data The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#'
#' @return A mully graph with the added drug layer
#' @export
#'
#' @examples
#' \dontrun{ 
#' data=readDBXML(DBXMLFilePath)
#' g=mully("DrugBank",direct=T)
#' g=addDBLayer(g,data,c("DB00001","DB06605"))
#' }
#' @import mully
#' @importFrom igraph V
addDBLayer<-function(g,data,drugList){
  dbmully=addLayer(g,"drugs")
  drugs=getDBDrug(data,drugList)
  interactions=getDBDrugInteractions(data,drugs$'drugbank_id')
  #Add Drugs' nodes
  message("Multipath: Adding Drugs Nodes")
  for (i in 1:dim(drugs)[1]) {
    progress(i,progress.bar = T)
    attrList=drugs[i,]
    attr=as.list(attrList)
    names(attr)=names(drugs)
    dbmully=mully::addNode(dbmully,nodeName = drugs$'drugbank_id'[i],layerName = "drugs",attributes = attr[-1])
  }
  
  message("Multipath: DONE - Drugs Nodes Added")
  
  #Add Drugs' interactions
  message("Multipath: Adding Drugs' Interactions")
  for (i in 1:dim(interactions)[1]) {
    progress(i,progress.bar = T)
    startName=V(dbmully)[which(V(dbmully)$name == interactions$'target_drugbank_id'[i])]$name
    endName=V(dbmully)[which(V(dbmully)$name == interactions$'drugbank_id'[i])]$name
    if(!is.null(startName) & !is.null(endName) & length(startName)!=0 & length(endName)!=0 )
      dbmully=mully::addEdge(dbmully,startName,endName,attributes = list(description=interactions$description[i]))
  }
  message("Multipath: DONE - Drugs' Interactions Added")
  return(dbmully)
}

#' Get the targets of given DrugBank drugs
#'
#' @param data  The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#'
#' @return A dataframe containing all information on the targets of the given drug list
#' @export
getDBTargets<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'drugbank_id'
  #Get Targets' list
  allTargets=dbparser:::targets()
  targets=allTargets[which(allTargets$'drugbank_id'%in%drugs),]
  
  #Get Targets' infos
  polypeptides=as.data.frame(data[["targets_polypeptides"]])
  infos=polypeptides[which(polypeptides$'target_id'%in%targets$'target_id'),]

  #Join the dataframes
  result=merge.data.frame(targets,infos)
  
  #Rename columns
  colnames(result)[which(colnames(result)=="drugbank_id")]<-"dbid"
  return(result)
}

#' Get the enzymes inhibited/induced or involved in metabolism by given DrugBank drugs
#'
#' @param data  The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#'
#' @return A dataframe containing all information on the enzymes inhibited/induced or involved in metabolism by the given drug list
#' @export
#'
getDBEnzymes<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'drugbank_id'
  #Get Enzymes' list
  allEnzymes=dbparser:::enzymes()
  enzymes=allEnzymes[which(allEnzymes$'drugbank_id'%in%drugs),]
  
  #Get Targets' infos
  polypeptides=dbparser:::enzymes_polypeptides()
  infos=polypeptides[which(polypeptides$'enzyme_id'%in%enzymes$'enzyme_id'),]

  #Join the dataframes
  result=merge.data.frame(enzymes,infos)
  
  #Rename columns
  colnames(result)[which(colnames(result)=="drugbank_id")]<-"dbid"
  return(result)
}


#' Get the transporter proteins involved in movement of given drugs across biological membranes
#'
#' @param data  The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#'
#' @return A dataframe containing all information on the transporter proteins involved in movement of the given drugs across biological membranes
#' @export
#'
getDBTransporters<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'drugbank_id'
  #Get Transporters' list
  allTransporters=dbparser:::transporters()
  transporters=allTransporters[which(allTransporters$'drugbank_id'%in%drugs),]
  
  #Get Transporters' infos
  polypeptides=dbparser:::transporters_polypeptides()
  infos=polypeptides[which(polypeptides$'transporter_id'%in%transporters$'transporter_id'),]
  
  #Join the dataframes
  result=merge.data.frame(transporters,infos)
  
  #Rename columns
  colnames(result)[which(colnames(result)=="drugbank_id")]<-"dbid"
  return(result)
}

#' Get the carrier proteins involved in movement of given drugs across biological membranes
#'
#' @param data  The dataframe containing the parsed information of DrugBank. This argument can be obtained using the function loadDBXML(DrugBankFile)
#' @param drugList The list of DrugBank Ids of the drugs. This argument can be either a string (one drug) or a list of strings (multiple drugs)
#'
#' @return A dataframe containing all information on the carrier proteins involved in movement of the given drugs across biological membranes
#' @export
#'
getDBCarriers<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'drugbank_id'
  #Get Carriers' list
  allCarriers=dbparser:::carriers()
  carriers=allCarriers[which(allCarriers$'drugbank_id'%in%drugs),]
  
  #Get Carriers' infos
  polypeptides=dbparser:::carriers_polypeptides()
  infos=polypeptides[which(polypeptides$'carrier_id'%in%carriers$'carrier_id'),]

  #Join the dataframes
  result=merge.data.frame(carriers,infos)
  
  #Rename columns
  colnames(result)[which(colnames(result)=="drugbank_id")]<-"dbid"
  return(result)
}