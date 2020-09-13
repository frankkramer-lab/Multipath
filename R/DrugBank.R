#' Load DrugBank XML file
#' 
#' @param file The path to the DrugBank XML file. This can be downloaded from the DrugBank official Website (drugbank.ca). An account with an institutional e-mail is required.
#'
#' @return A dataframe containing the parsed information from DrugBank. This can be used to extract any additional information on the DrugBank entries
#' @export
#'
#' @note
#' This function should be called before using any function to query the DrugBank database. Since the parsing of DrugBank takes time, this function should only be called once.
#'
loadDBXML<-function(file){
  read_drugbank_xml_db(file)
  data=run_all_parsers()
  return(data)
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
#' data=loadDBXML(DrugBankFilePath)
#' getDBDrug(data, "DB00001")
getDBDrug<-function(data,drug){
  drugs=as.data.frame(data[[1]])
  entry=drugs[which(drugs$'primary_key'%in%drug),]
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
#' getDBDrugInteractions(data,"DB06605")
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
  interactions=as.data.frame(data[["interactions"]])
  drugInteractions=interactions[which(interactions$'drugbank-id'%in%drugs & interactions$'parent_key'%in%drugs),]
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
#' addDBLayer(mully("DrugBank",direct=T),data,c("DB00001","DB06605"))
addDBLayer<-function(g,data,drugList){
  dbmully=addLayer(g,"drugs")
  drugs=getDBDrug(data,drugList)
  interactions=getDBDrugInteractions(data,drugs$'primary_key')
  #Add Drugs' nodes
  message("Multipath: Adding Drugs Nodes")
  for (i in 1:dim(drugs)[1]) {
    progress(i,progress.bar = T)
    attrList=drugs[i,]
    attr=as.list(attrList)
    names(attr)=names(drugs)
    dbmully=mully::addNode(dbmully,nodeName = drugs$'primary_key'[i],layerName = "drugs",attributes = attr[-1])
  }
  
  message("Multipath: DONE - Drugs Nodes Added")
  
  #Add Drugs' interactions
  message("Multipath: Adding Drugs' Interactions")
  for (i in 1:dim(interactions)[1]) {
    progress(i,progress.bar = T)
    startName=V(dbmully)[which(V(dbmully)$name == interactions$'drugbank-id'[i])]$name
    endName=V(dbmully)[which(V(dbmully)$name == interactions$'parent_key'[i])]$name
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
#' @return A dataframe containing all information on the tragets of the given drug list
#' @export
#'
#' @examples
getDBTargets<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'primary_key'
  #Get Targets' list
  allTargets=as.data.frame(data[["targets"]])
  targets=allTargets[which(allTargets$'parent_key'%in%drugs),]
  
  #Get Targets' infos
  polypeptides=as.data.frame(data[["targets_polypeptides"]])
  infos=polypeptides[which(polypeptides$'parent_id'%in%targets$id),]

  #Join the dataframes
  result=merge.data.frame(targets,infos,by.x=c("id","name","organism"),by.y=c("parent_id","name","organism"))
  
  #Rename columns
  colnames(result)[which(colnames(result)=="parent_key")]<-"dbid"
  colnames(result)[which(colnames(result)=="id.y")]<-"upid"
  colnames(result)[which(colnames(result)=="id")]<-"dbproteinid"
  colnames(result)[which(colnames(result)=="name")]<-"name"
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
#' @examples
getDBEnzymes<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'primary_key'
  #Get Enzymes' list
  allEnzymes=as.data.frame(data[["enzymes"]])
  enzymes=allEnzymes[which(allEnzymes$'parent_key'%in%drugs),]
  
  #Get Targets' infos
  polypeptides=as.data.frame(data[["enzymes_polypeptides"]])
  infos=polypeptides[which(polypeptides$'parent_id'%in%enzymes$id),]

  #Join the dataframes
  result=merge.data.frame(enzymes,infos,by.x=c("id","name","organism"),by.y=c("parent_id","name","organism"))
  
  #Rename columns
  colnames(result)[which(colnames(result)=="parent_key")]<-"dbid"
  colnames(result)[which(colnames(result)=="id.y")]<-"upid"
  colnames(result)[which(colnames(result)=="id")]<-"dbproteinid"
  colnames(result)[which(colnames(result)=="name")]<-"name"
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
#' @examples
getDBTransporters<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'primary_key'
  #Get Transporters' list
  allTransporters=as.data.frame(data[["transporters"]])
  transporters=allTransporters[which(allTransporters$'parent_key'%in%drugs),]
  
  #Get Transporters' infos
  polypeptides=as.data.frame(data[["transporters_polypeptides"]])
  infos=polypeptides[which(polypeptides$'parent_id'%in%transporters$id),]
  
  #Join the dataframes
  result=merge.data.frame(transporters,infos,by.x=c("id","name","organism"),by.y=c("parent_id","name","organism"))
  
  #Rename columns
  colnames(result)[which(colnames(result)=="parent_key")]<-"dbid"
  colnames(result)[which(colnames(result)=="id.y")]<-"upid"
  colnames(result)[which(colnames(result)=="id")]<-"dbproteinid"
  colnames(result)[which(colnames(result)=="name")]<-"name"
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
#' @examples
getDBCarriers<-function(data,drugList){
  drugs=getDBDrug(data,drugList)$'primary_key'
  #Get Carriers' list
  allCarriers=as.data.frame(data[["carriers"]])
  carriers=allCarriers[which(allCarriers$'parent_key'%in%drugs),]
  
  #Get Carriers' infos
  polypeptides=as.data.frame(data[["carriers_polypeptides"]])
  infos=polypeptides[which(polypeptides$'parent_id'%in%carriers$id),]

  #Join the dataframes
  result=merge.data.frame(carriers,infos,by.x=c("id","name","organism"),by.y=c("parent_id","name","organism"))
  
  #Rename columns
  colnames(result)[which(colnames(result)=="parent_key")]<-"dbid"
  colnames(result)[which(colnames(result)=="id.y")]<-"upid"
  colnames(result)[which(colnames(result)=="id")]<-"dbproteinid"
  colnames(result)[which(colnames(result)=="name")]<-"name"
  return(result)
}

# addDBTargetLayer<-function(g,data,drugList){
#   dbmully=addLayer(g,"targets")
#   drugs=getDBDrug(data,drugList)
# 
#   targets=data[["drug_targets"]]
#   #Add Targets' nodes and relations to Drugs
#   for (i in 1:dim(targets)[1]) {
#     targetName=targets$name[i]
#     idTarget=getIDNode(dbmully,nameNode = targetName)
#     if(is.null(idTarget)){
#       attr=list(beid=targets$id[i],organism=targets$organism[i])
#       dbmully=addNode(dbmully,nodeName = targets$name[i],layerName = "targets",attributes = attr)
#       idTarget=getIDNode(dbmully,nameNode = targetName)
#     }
#     nameDrug=V(dbmully)[which(V(dbmully)$dbid == targets$parent_key[i])]$name
#     if(!is.null(nameDrug)){
#       dbmully=addEdge(dbmully,nameDrug,targetName,attributes = list(type="bind"))
#     }
#   }
#   return(dbmully)
# }

# dbmully=mully("DrugBank",direct = T)
# dbmully=addLayer(dbmully,c("drugs","targets"))
# drugs=data[[1]]
# 
# drugBank=data
# for(i in 1:72){
#   any(duplicated(data[[i]]$`drugbank-id`))
# }