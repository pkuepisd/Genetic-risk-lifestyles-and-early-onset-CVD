# Incidence rate change
library(data.table)
library(survival)
library(foreach)
library(doParallel)
library(lmtest)
library(boot)

Args<- commandArgs()
ep<- Args[3]

setwd("/public/home/sundong/PhD/Project3")
load(paste0("GxE_",ep,"_AM4.RData"))

group<- "GR_LF7"
adj<- "+national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array+education_3g2+education_3g3+marital_2g+is_female"

# boot
IR_boot<- function(data,inds) {
	DF<- as.data.frame(data)[inds,]
	DF$group<- DF[,group]

	DF<- within(DF,{
			ep_early<- NA
			ep_early[is_female==0 & age_ep<55 & ep==1]<- 1
			ep_early[is_female==1 & age_ep<65 & ep==1]<- 1
	})
	DF$ep_early[is.na(DF$ep_early)]<- 0

	DF<- within(DF,{
		ep_late<- NA
		ep_late[is_female==0 & age_ep>=55 & ep==1]<- 1
		ep_late[is_female==1 & age_ep>=65 & ep==1]<- 1
	})
	DF$ep_late[is.na(DF$ep_late)]<- 0

	levels<- length(table(DF$group))
	level_label<- attributes(DF$group)$levels

	early_onset_data<- subset(DF,(is_female==0 & age_enter<55) | (is_female==1 & age_enter<65))
	early_onset_data<- within(early_onset_data,{
		age_ep_early<- age_ep
		age_ep_early[is_female==0 & age_ep>=55]<- 55
		age_ep_early[is_female==1 & age_ep>=65]<- 65
		})
	early_onset_data$fp<- with(early_onset_data,age_ep_early-age_enter)

	formula<- paste0("Surv(fp,ep_early)~","group+age_enter",adj)
	fit<- coxph(as.formula(formula),data=early_onset_data)

	newdf<- with(early_onset_data,data.frame(fp=max(fp),age_enter=mean(age_enter),ep_early=1,ep_late=1,group=level_label))
	newdf2<- t(mapply(mean,subset(DF,select=-c(age_enter,ep_early,ep_late,group))))
	newdf<- cbind(newdf,newdf2)
	pre<- predict(fit,newdf,type="expected")
	result<- pre/max(early_onset_data$fp)*100000
	final<- result

	late_onset_data<- subset(DF,(is_female==0 & age_ep>=55) | (is_female==1 & age_ep>=65))
	late_onset_data<- within(late_onset_data,{
		age_enter_late<- age_enter
		age_enter_late[is_female==0 & age_enter<55]<- 55
		age_enter_late[is_female==1 & age_enter<65]<- 65
	})
	late_onset_data$fp<- with(late_onset_data,age_ep-age_enter_late)
	formula<- paste0("Surv(fp,ep_late)~","group+age_enter_late",adj)
	fit<- coxph(as.formula(formula),data=late_onset_data)

	newdf<- with(late_onset_data,data.frame(fp=max(fp),age_enter_late=mean(age_enter_late),ep_early=1,ep_late=1,group=level_label))
	newdf2<- t(mapply(mean,subset(DF,select=-c(ep_early,ep_late,group))))
	newdf<- cbind(newdf,newdf2)

	pre<- predict(fit,newdf,type="expected")
	result<- pre/max(late_onset_data$fp)*100000
	final<- c(final,result)
	return(final)
}

crude<- boot(full_prs_testing,IR_boot,R=999)
result<- c()
for (j in 1:32) {
	est_j<- sprintf("%0.2f",crude$t0[j])
	lci_j<- sprintf("%0.2f",sort(crude$t[,j])[25])
	uci_j<- sprintf("%0.2f",sort(crude$t[,j])[975])
	est_ci<- paste0(est_j," (",lci_j,"-",uci_j,")")
	result<- c(result,est_ci)
}