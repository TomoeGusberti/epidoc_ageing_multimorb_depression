parTableToCSV<-function(lavObj,path){
  
  if(inherits(lavObj,'lavaan.mi')){
    # check if group analysis
    listcall<-lavObj@lavListCall
    group<-("group"%in%names(listcall))
    if(group){
      # group
      labels<-lavObj@Data@group.label
      labelsDf<-data.frame(
        group=c(1:length(labels)),
        Group=labels
      )
    }
    partable<-parameterEstimates.mi(lavObj,standardized =TRUE)
  }else{
    listcall<-lavObj@call
    group<-("group"%in%names(listcall))
    if(group){
      # group
      labels<-lavObj@Data@group.label
    }
    partable<-parameterEstimates(lavObj,standardized =TRUE)
  }
  # save parameters
  if(group){
    labelsDf<-data.frame(
      group=c(1:length(labels)),
      Group=labels
    )
    partable0<-partable%>%mutate(
      label2=paste0(lhs,op,rhs),
      estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
      CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
    )
    partable0=left_join(partable0,labelsDf)
    partable<-reshape(partable0%>%filter(!block==0), idvar = c('label2'),
                       timevar = "Group", direction = "wide",sep='_')
    gPartable<-partable%>%
      select(label2,
             names(partable)[sapply(names(partable),function(x){grepl('estimate',x)})],
             names(partable)[sapply(names(partable),function(x){grepl('CI',x)})]
      )
    
    
    write.csv(gPartable,
              file=paste0(path,'/parTable_Mgroup_',listcall$group,'.csv')
    )
    
    if ('label'%in%names(partable0)){
      partable_specific<-partable0%>%filter(block==0)%>%mutate(Group2=as.numeric(substr(label,2,2)),
                                                               label=substr(label,3,nchar(label)))%>%
        select(Group2,label,estimate_p, CI)%>%
        left_join(labelsDf,by=c('Group2'='group'))
      print(partable_specific)
      partable_specific<-reshape(partable_specific%>%select(-Group2), idvar = c('label'),
                                 timevar = "Group", direction = "wide",sep='_')
      write.csv(partable_specific,
                file=paste0(path,'/parTable_Mgroup_',listcall$group,'_specifics.csv'))
      
      partable_spec2=partable0%>%filter(op==':=')%>%select(label2,estimate_p, CI)
      write.csv(partable_spec2,
                file=paste0(path,'/parTable_Mgroup_',listcall$group,'_specifiedVars.csv'))
    }
    
    res<-list(partable=partable0,gPartable=gPartable)
  } else{
    # no group
    partable<-partable%>%mutate(
      label2=paste0(lhs,op,rhs),
      estimate_p=paste0(round(est,3),' (',round(pvalue,3),')'),
      CI=paste0('[',round(ci.lower,3),':',round(ci.upper,3),']')
    )%>%select(-c(lhs,op,rhs))
    write.csv(partable,file=paste0(path,'/parTable_avgModel.csv'))
    res<-list(partable=partable)
  }
  return(res)
  
  
  
}