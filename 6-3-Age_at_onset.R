# Associations of genetic risk and lifestyles with different age-onset of CVD (R file) 
library(data.table)
library(survival)
library(foreach)
library(doParallel)
library(lmtest)

cl<- makeCluster(20)
registerDoParallel(cl)
 
Args<- commandArgs()
ep<- Args[3]
index<- as.numeric(Args[4])

setwd("/public/home/sundong/PhD/Project3")
load(paste0("GxE_",ep,"_AM4.RData"))

GRadj<- "national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array"
baseadj<- "education_3g+marital_2g+is_female"
region<- "region_code"
varlist<- c("MetaPRS","GR","LF7","GR_LF7","GRT","GRT_LF7","cLF7")

adjlist<- c(
	rep(paste0(GRadj,"+",baseadj),2),
	paste0(region,"+",baseadj),
	rep(paste0(GRadj,"+",baseadj),3),
	paste0(region,"+",baseadj)
)

group<- varlist[index]
adj<- adjlist[index]

TD_effect<- function(group,adj,cutpoints,data,sex="All"){
	DF<- as.data.frame(data)
	if (sex=="Men") DF<- DF[DF$is_female==0,]
	if (sex=="Women") DF<- DF[DF$is_female==1,]

	DF$group<- DF[,group]

	formula1<- paste0("Surv(age_enter,age_ep,ep)~","group+",adj)
	PH_table_before<- cox.zph(coxph(as.formula(formula1),data=DF))$table
	PH_test_before<- sprintf("%0.3f",PH_table_before[1,3])
	if (PH_test_before=="0.000") PH_test_before<- "<0.001"

	if (is.factor(DF$group)) {
		DF$group<- factor(DF$group,levels=rev(attributes(DF$group)$levels))
		levels<- length(table(DF$group))
		level_label<- attributes(DF$group)$levels
	} 

	split_data<- survSplit(Surv(age_enter,age_ep,ep)~.,data=DF,cut=cutpoints,episode="tgroup")
	formula2<- paste0("Surv(age_enter,age_ep,ep)~","group:strata(tgroup)+",adj)
	tdcox<- coxph(as.formula(formula2),data=split_data)

	PH_table_after<- cox.zph(tdcox)$table
	PH_test_after<- sprintf("%0.3f",PH_table_after[(nrow(PH_table_after)-1),3])
	if (PH_test_after==0) PH_test_after<- "<0.001"

	curve<- survfit(Surv(age_enter,age_ep,ep)~ 1,data=split_data)
	event_list<- summary(curve,time=c(cutpoints,max(DF$age_ep)))
	fp<- split_data$age_ep- split_data$age_enter

	N_cut<- length(cutpoints)+1
	for (i in 1:N_cut) {
		if (!is.factor(DF$group)) {
			hr<-sprintf("%0.2f",exp(rev(coefficients(tdcox))[N_cut-i+1]))
	        lci<-sprintf("%0.2f",exp(rev(confint(tdcox)[,1])[N_cut-i+1]))
	        uci<-sprintf("%0.2f",exp(rev(confint(tdcox)[,2])[N_cut-i+1]))
	        v<- paste0(hr," ","(",lci,"-",uci,")")
		} else {
	        v<- c("Reference")
	        for (j in 2:levels) {
	        	hr<-sprintf("%0.2f",exp(rev(coefficients(tdcox))[(N_cut-i)*levels+j]))
		        lci<-sprintf("%0.2f",exp(rev(confint(tdcox)[,1])[(N_cut-i)*levels+j]))
		        uci<-sprintf("%0.2f",exp(rev(confint(tdcox)[,2])[(N_cut-i)*levels+j]))
		        hr_ci<- paste0(hr," ","(",lci,"-",uci,")")
		        v<- c(v,hr_ci)
	        }
	    }
        case<- event_list$n.event[i]
        rate<- sprintf("%0.2f",case/sum(fp[split_data$tgroup==i])*1000)
        case<- format(case,big.mark=",")
        if (i==1) label<- paste0("-",cutpoints[1])
        if (i!=1 & i!=(length(cutpoints)+1)) label<- paste0(cutpoints[i-1],"-",cutpoints[i])
        if (i==(length(cutpoints)+1)) label<- paste0(cutpoints[length(cutpoints)],"-")

        tgroup_result<- c(label,case,rate,v,PH_test_before,PH_test_after)
        if (i==1) final<- tgroup_result
        if (i!=1) final<- rbind(final,tgroup_result)
	}
	final<- rbind(final," ")
	return(final)
}

cutpoints<- seq(60,80,10)
TD_result<- TD_effect(group,adj,cutpoints,full_prs_testing)
TD_result_male<- TD_effect(group,adj,cutpoints,full_prs_testing,"Men")
TD_result_female<- TD_effect(group,adj,cutpoints,full_prs_testing,"Women")

SP_effect<- function(group,adj,age_cut_men=55,age_cut_women=65,data){
	DF<- as.data.frame(data)
	DF$group<- DF[,group]

	DF<- within(DF,{
		ep_early<- NA
		ep_early[is_female==0 & age_ep<age_cut_men & ep==1]<- 1
		ep_early[is_female==1 & age_ep<age_cut_women & ep==1]<- 1
	})
	DF$ep_early[is.na(DF$ep_early)]<- 0

	DF<- within(DF,{
		ep_late<- NA
		ep_late[is_female==0 & age_ep>=age_cut_men & ep==1]<- 1
		ep_late[is_female==1 & age_ep>=age_cut_women & ep==1]<- 1
	})
	DF$ep_late[is.na(DF$ep_late)]<- 0

	if (is.factor(DF$group)) {
		levels<- length(table(DF$group))
		level_label<- attributes(DF$group)$levels
	} 

	for (time in c("early","late")){
		formula<- paste0("Surv(age_enter,age_ep,ep_",time,")~","group+",adj)
		fit<- coxph(as.formula(formula),data=DF)

		if (!is.factor(DF$group)) {
			hr<-sprintf("%0.2f",exp(coefficients(fit)[1]))
	        lci<-sprintf("%0.2f",exp(confint(fit)[1,1]))
	        uci<-sprintf("%0.2f",exp(confint(fit)[1,2]))
	        v<- paste0(hr," ","(",lci,"-",uci,")")
		} else {
	        v<- c("Reference")
	        for (j in 1:(levels-1)) {
	        	hr<-sprintf("%0.2f",exp(coefficients(fit)[j]))
		        lci<-sprintf("%0.2f",exp(confint(fit)[j,1]))
		        uci<-sprintf("%0.2f",exp(confint(fit)[j,2]))
		        hr_ci<- paste0(hr," ","(",lci,"-",uci,")")
		        v<- c(v,hr_ci)
	        }
	    }

	    if (group=="GR_LF7") {
	    	formula1<- paste0("Surv(age_enter,age_ep,ep_",time,")~","GR+LF7+",adj)
	    	fit1<- coxph(as.formula(formula1),data=DF)
	    	LR_test<- sprintf("%0.3f",lrtest(fit1,fit)[2,5])
			if (LR_test=="0.000") LR_test<- "<0.001"
	    }
	    if (group=="GRT_LF7") {
	    	formula1<- paste0("Surv(age_enter,age_ep,ep_",time,")~","GRT+LF7+",adj)
	    	fit1<- coxph(as.formula(formula1),data=DF)
	    	LR_test<- sprintf("%0.3f",lrtest(fit1,fit)[2,5])
			if (LR_test=="0.000") LR_test<- "<0.001"
	    }

		event<- nrow(DF[DF[,paste0("ep_",time)]==1,])
		if (time=="early") {
			early_onset_data<- subset(DF,(is_female==0 & age_enter<age_cut_men) | (is_female==1 & age_enter<age_cut_women))
			early_onset_data<- within(early_onset_data,{
				age_ep_early<- age_ep
				age_ep_early[is_female==0 & age_ep>=55]<- 55
				age_ep_early[is_female==1 & age_ep>=65]<- 65
				})
			fp<- with(early_onset_data,age_ep_early-age_enter)
			label<- paste0("Male<",age_cut_men," | Female<",age_cut_women)
		} else {
			fp<- with(DF,age_ep-age_enter)
			label<- paste0("Male>=",age_cut_men," | Female>=",age_cut_women)
		}

		rate<- sprintf("%0.2f",event/sum(fp)*1000)
		event<- format(event,big.mark=",")
		PH_test<- sprintf("%0.3f",cox.zph(fit)$table[1,3])

        tgroup_result<- c(label,event,rate,v,PH_test," ")
        if (time=="early") final<- tgroup_result
        if (time!="early") final<- rbind(final,tgroup_result)
	}
	final<- rbind(final," ")
	return(final)
}

SP_result<- SP_effect(group,adj,55,65,full_prs_testing)