ptrn="^Breast"


getCasesIDList<-function(pattern){
  caseidlist=get_cancer_studies()
  caselist=caseidlist[str_detect(caseidlist$name,pattern),]$cancer_study_id
  return(caselist)
}

getTCGAData<-function(idList){
  tcgadata=get_case_lists(idList[1])
  rows=length(idlist)
  for (i in 2:rows) {
    tcgadata=rbind.data.frame(tcgadata,get_case_lists(idlist[i]))
  }
  return(tcgadata)
}


# allbreastalp20=unlist(strsplit(cancerdata[1,5]," "))
