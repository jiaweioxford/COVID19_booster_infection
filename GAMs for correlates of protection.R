####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# for analyses of correlates of protection against SARS-Cov-2 infection.                           #
# Accompanying paper: Protection against SARS-CoV-2 Omicron BA.4/5 variant following booster       #
# vaccination or breakthrough infection in the UK                                                  #
####################################################################################################

library(mgcv)
library(dplyr)

## model fitting
m <- bam(result_mk ~ s(assay_21_59,k=200)+s(dur_vac_inf,k=20)+
           s(age_at_visit,k=20)+s(study_day,k=20)+
           ti(assay_21_59,dur_vac_inf,k=c(200,20)) +
           ti(assay_21_59,age_at_visit,k=c(200,20)) + variant+
           ti(study_day, age_at_visit, k=c(20,20), by=region_model) + 
           ethnicity_wo + sex + region_model + 
           imd_pc_new + multigen + hhsizegroup + rural_urban_class +
           ever_care_home_worker + patient_facing_clean_ever + 
           ever_personfacing_socialcare + ever_lthc + health_condition+
           visit_freq + smoke_now + contact_hospital + contact_carehome, data = data, 
         family = binomial,method="fREML",discrete=TRUE,nthreads=8)


##function to extract odds ratio
getor_vaccine_type <- function(newdata = newdata, contrast = c(2,8), model=model) {
  newdata$assay_21_59 <- c(contrast[1],contrast[2])
  Xp <- predict(model, newdata = newdata, type="lpmatrix")
  diff <- (Xp[2,]-Xp[1,])
  dly <- t(diff)%*%coef(model)
  se.dly <- sqrt(t(diff)%*%vcov(model)%*%diff)
  c(exp(dly), exp(dly -2*se.dly), exp(dly + 2*se.dly))
}




## prediction for vaccination group

age_list=seq(20,80,5)
results_list=list()
j=1
for(age in age_list){
  newdata = data.frame(age_at_visit=age, sex=2,ethnicity_wo=1,study_day=400, region_model=8,imd_pc_new=60,multigen=1,
                       hhsizegroup=2,rural_urban_class=2, ever_care_home_worker=1,patient_facing_clean_ever=1, ever_personfacing_socialcare=1,
                       ever_lthc=1,visit_freq=1,smoke_now=1,contact_hospital=1,contact_carehome=1,dur_vac_inf=210,variant="Vaccination",health_condition="No")
  newdata <- rbind(newdata,newdata)
  
  dat_ors_vac <- data.frame(assay_round = c(16,seq(20,8000,20)))
  dat_ors_vac$or <- NA
  dat_ors_vac$ll <- NA
  dat_ors_vac$ul <- NA
  
  
  for (i in 1:nrow(dat_ors_vac)) {
    temp <- getor_vaccine_type(newdata=newdata,  contrast=c(16,dat_ors_vac$assay_round[i]),
                               model = m)
    dat_ors_vac$or[i] <- temp[1]
    dat_ors_vac$ll[i] <- temp[2]
    dat_ors_vac$ul[i] <- temp[3]
  }
  
  dat_ors_vac=dat_ors_vac%>%mutate(VE=100*(1-or),VE_ll=100*(1-ul),VE_ul=100*(1-ll))
  dat_ors_vac$variant="Vaccination"
  dat_ors_vac$age=age
  
  results_list[[j]]<-dat_ors_vac
  j=j+1
}

results_vac=do.call(rbind,results_list)



## prediction for Pre-Alpha/Alpha group

age_list=seq(20,80,5)
results_list=list()
j=1
for(age in age_list){
  newdata = data.frame(age_at_visit=age, sex=2,ethnicity_wo=1,study_day=400, region_model=8,imd_pc_new=60,multigen=1,
                       hhsizegroup=2,rural_urban_class=2, ever_care_home_worker=1,patient_facing_clean_ever=1, ever_personfacing_socialcare=1,
                       ever_lthc=1,visit_freq=1,smoke_now=1,contact_hospital=1,contact_carehome=1,dur_vac_inf=210,health_condition="No")
  newdata <- rbind(newdata,newdata)
  newdata$variant <- c("Vaccination","Pre-Alpha/Alpha")
  
  dat_ors_vac <- data.frame(assay_round = seq(80,8000,20))
  dat_ors_vac$or <- NA
  dat_ors_vac$ll <- NA
  dat_ors_vac$ul <- NA
  
  
  for (i in 1:nrow(dat_ors_vac)) {
    temp <- getor_vaccine_type(newdata=newdata,  contrast=c(16,dat_ors_vac$assay_round[i]),
                               model = m)
    dat_ors_vac$or[i] <- temp[1]
    dat_ors_vac$ll[i] <- temp[2]
    dat_ors_vac$ul[i] <- temp[3]
  }
  
  dat_ors_vac=dat_ors_vac%>%mutate(VE=100*(1-or),VE_ll=100*(1-ul),VE_ul=100*(1-ll))
  dat_ors_vac$variant="Pre-alpha/Alpha"
  dat_ors_vac$age=age
  
  results_list[[j]]<-dat_ors_vac
  j=j+1
}

results_alpha=do.call(rbind,results_list)




## prediction for Delta group

age_list=seq(20,80,5)
results_list=list()
j=1
for(age in age_list){
  newdata = data.frame(age_at_visit=age, sex=2,ethnicity_wo=1,study_day=400, region_model=8,imd_pc_new=60,multigen=1,
                       hhsizegroup=2,rural_urban_class=2, ever_care_home_worker=1,patient_facing_clean_ever=1, ever_personfacing_socialcare=1,
                       ever_lthc=1,visit_freq=1,smoke_now=1,contact_hospital=1,contact_carehome=1,dur_vac_inf=180,health_condition="No")
  newdata <- rbind(newdata,newdata)
  newdata$variant <- c("Vaccination","Delta")
  
  dat_ors_vac <- data.frame(assay_round = seq(100,8000,20))
  dat_ors_vac$or <- NA
  dat_ors_vac$ll <- NA
  dat_ors_vac$ul <- NA
  
  
  for (i in 1:nrow(dat_ors_vac)) {
    temp <- getor_vaccine_type(newdata=newdata,  contrast=c(16,dat_ors_vac$assay_round[i]),
                               model = m)
    dat_ors_vac$or[i] <- temp[1]
    dat_ors_vac$ll[i] <- temp[2]
    dat_ors_vac$ul[i] <- temp[3]
  }
  
  dat_ors_vac=dat_ors_vac%>%mutate(VE=100*(1-or),VE_ll=100*(1-ul),VE_ul=100*(1-ll))
  dat_ors_vac$variant="Delta"
  dat_ors_vac$age=age
  
  results_list[[j]]<-dat_ors_vac
  j=j+1
}

results_delta=do.call(rbind,results_list)



## prediction for BA.1 group

age_list=seq(20,80,5)
results_list=list()
j=1
for(age in age_list){
  newdata = data.frame(age_at_visit=age, sex=2,ethnicity_wo=1,study_day=400, region_model=8,imd_pc_new=60,multigen=1,
                       hhsizegroup=2,rural_urban_class=2, ever_care_home_worker=1,patient_facing_clean_ever=1, ever_personfacing_socialcare=1,
                       ever_lthc=1,visit_freq=1,smoke_now=1,contact_hospital=1,contact_carehome=1,dur_vac_inf=180,health_condition="No")
  newdata <- rbind(newdata,newdata)
  newdata$variant <- c("Vaccination","BA.1")
  
  dat_ors_vac <- data.frame(assay_round = seq(140,8000,20))
  dat_ors_vac$or <- NA
  dat_ors_vac$ll <- NA
  dat_ors_vac$ul <- NA
  
  
  for (i in 1:nrow(dat_ors_vac)) {
    temp <- getor_vaccine_type(newdata=newdata,  contrast=c(16,dat_ors_vac$assay_round[i]),
                               model = m)
    dat_ors_vac$or[i] <- temp[1]
    dat_ors_vac$ll[i] <- temp[2]
    dat_ors_vac$ul[i] <- temp[3]
  }
  
  dat_ors_vac=dat_ors_vac%>%mutate(VE=100*(1-or),VE_ll=100*(1-ul),VE_ul=100*(1-ll))
  dat_ors_vac$variant="BA.1"
  dat_ors_vac$age=age
  
  results_list[[j]]<-dat_ors_vac
  j=j+1
}

results_ba1=do.call(rbind,results_list)

## prediction for BA.2 group

age_list=seq(20,80,5)
results_list=list()
j=1
for(age in age_list){
  newdata = data.frame(age_at_visit=age, sex=2,ethnicity_wo=1,study_day=400, region_model=8,imd_pc_new=60,multigen=1,
                       hhsizegroup=2,rural_urban_class=2, ever_care_home_worker=1,patient_facing_clean_ever=1, ever_personfacing_socialcare=1,
                       ever_lthc=1,visit_freq=1,smoke_now=1,contact_hospital=1,contact_carehome=1,dur_vac_inf=120,health_condition="No")
  newdata <- rbind(newdata,newdata)
  newdata$variant <- c("Vaccination","BA.2")
  
  dat_ors_vac <- data.frame(assay_round = seq(200,8000,20))
  dat_ors_vac$or <- NA
  dat_ors_vac$ll <- NA
  dat_ors_vac$ul <- NA
  
  
  for (i in 1:nrow(dat_ors_vac)) {
    temp <- getor_vaccine_type(newdata=newdata,  contrast=c(16,dat_ors_vac$assay_round[i]),
                               model = m)
    dat_ors_vac$or[i] <- temp[1]
    dat_ors_vac$ll[i] <- temp[2]
    dat_ors_vac$ul[i] <- temp[3]
  }
  
  dat_ors_vac=dat_ors_vac%>%mutate(VE=100*(1-or),VE_ll=100*(1-ul),VE_ul=100*(1-ll))
  dat_ors_vac$variant="BA.2"
  dat_ors_vac$age=age
  
  results_list[[j]]<-dat_ors_vac
  j=j+1
}

results_ba2=do.call(rbind,results_list)

