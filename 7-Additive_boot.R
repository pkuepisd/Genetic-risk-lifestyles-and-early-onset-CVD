# Additive interaction
library(boot)
library(data.table)
library(survival)
library(openxlsx)

for (ep in c("cad","is","hs")) {
	setwd("/public/home/sundong/PhD/Project3")
	load(paste0("GxE_",ep,"_AM4.RData"))
	GRadj<- "national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array"
	baseadj<- "education_3g+marital_2g+is_female"
	region<- "region_code"

	group<- "GR_LF7"
	adj<- paste0(GRadj,"+",baseadj)

	DF<- as.data.frame(full_prs_testing)
	DF$group<- DF[,group]
	age_cut_men<- 55
	age_cut_women<- 65

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
	levels<- length(table(DF$group))
	level_label<- attributes(DF$group)$levels

	Additive_boot<- function(data,inds) {
		DF<- data[inds,]
		for (time in c("early","late")){
			formula<- paste0("Surv(age_enter,age_ep,ep_",time,")~","group+",adj)
			fit<- coxph(as.formula(formula),data=DF)

			# RERI
			RERI<- exp(coef(fit)[8])-exp(coef(fit)[2])-exp(coef(fit)[6])+1
			# AP
			AP<-  (RERI/exp(coef(fit)[8]))*100
			# AP*
			# G only
			G<- ((exp(coef(fit)[6])-1)/(exp(coef(fit)[8])-1))*100			
			# LF only
			LF<- ((exp(coef(fit)[2])-1)/(exp(coef(fit)[8])-1))*100		
			# Interaction
			Int<- (RERI/(exp(coef(fit)[8])-1))*100
			result_high<- c(RERI,AP,G,LF,Int)

			# RERI
			RERI_HI<- exp(coef(fit)[7])-exp(coef(fit)[1])-exp(coef(fit)[6])+1

			# RERI
			RERI_IU<- exp(coef(fit)[5])-exp(coef(fit)[3])-exp(coef(fit)[2])+1

			# RERI
			RERI_II<- exp(coef(fit)[4])-exp(coef(fit)[3])-exp(coef(fit)[1])+1

			result<- c(result_high,RERI_HI,RERI_IU,RERI_II)
			if (time=="early") final<- result else final<- c(final,result)
		}
		return(final)
	}
	set.seed(1234)
	crude<- boot(DF,Additive_boot,R=999)
	result<- c()
	for (j in 1:16) {
		if (j>=2 & j<=5) {
			est_j<- sprintf("%0.1f",crude$t0[j])
			lci_j<- sprintf("%0.1f",sort(crude$t[,j])[25])
			uci_j<- sprintf("%0.1f",sort(crude$t[,j])[975])
		} else {
			est_j<- sprintf("%0.2f",crude$t0[j])
			lci_j<- sprintf("%0.2f",sort(crude$t[,j])[25])
			uci_j<- sprintf("%0.2f",sort(crude$t[,j])[975])
		}
		est_ci<- paste0(est_j," (",lci_j,"-",uci_j,")")
		result<- c(result,est_ci)
	}
	if (ep=="cad") final<- result else final<- rbind(final,result)
}