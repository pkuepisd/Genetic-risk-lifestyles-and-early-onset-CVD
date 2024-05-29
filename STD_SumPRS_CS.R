library(data.table)
library(plyr)
library(doParallel)
library(foreach)
library(survival)
cl<- makeCluster(20)
registerDoParallel(cl)
Args<-commandArgs()

gwas<-Args[3]

prspath<- paste0("~/PhD/Project1/GWASdata/",gwas,"/",gwas,".org.stdqc1_chr/CS")

for (phi in c("auto","1","1e-2","1e-4","1e-6")) {
	chr1 <- fread(paste0(prspath,"/chr1_",phi,".prs.profile"),header=T)[,c(1,6)]
	chrname<- paste0("prs",1:22)
	PRS<- chr1
	colnames(PRS)[1:2]<-c("id",chrname[1])
	if (phi=="auto") {
		count<- nrow(fread(paste0(prspath,"/chr1_",phi,".txt"),header=FALSE))
	}
	for (i in 2:22){
		df<- fread(paste0(prspath,"/chr",i,"_",phi,".prs.profile"),header=T)[,c(1,6)]
		colnames(df)[1:2]<- c("id",chrname[i])
		PRS<- merge(PRS,df,by="id")
		if (phi=="auto") {
			count_i<- nrow(fread(paste0(prspath,"/chr",i,"_",phi,".txt"),header=FALSE))
			count<- c(count,count_i)
		}
	}
	if (phi=="auto") {
		count<- sum(count)
	}
	PRS[,phi] <- rowSums(PRS[,2:23])
	PRS<- subset(PRS,select=c("id",phi))
	if (phi=="auto") {
		ALL<- PRS
	} else {
		ALL<- merge(ALL,PRS,by="id",all.x=TRUE)
	}
}
colnames(ALL)<- c("id","phi_auto","phi_1","phi_2","phi_3","phi_4")
fwrite(ALL,paste0(prspath,"/id_CSprs_ALL.txt"),sep=" ", na="NA",col.names=TRUE)

load(prspath,"/",gwas,"_full.training.RData")
load(prspath,"/",gwas,"_full.testing.RData")

eplist<- c("cad","is","hs")
for (ep in eplist){
	cox_result<- function(ep,ex,trait,count,df){
		if (ex=="phi_auto") {
			phi<- "auto"
		}
		if (ex=="phi_1") {
			phi<- "1"
		}
		if (ex=="phi_2") {
			phi<- "1e-2"
		}
		if (ex=="phi_3") {
			phi<- "1e-4"
		}
		if (ex=="phi_4") {
			phi<- "1e-6"
		}
	    formula<- paste0("Surv(age_enter,age_",ep,",",ep,")~",ex,"+national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array+strata(is_female)")
	    a<-coxph(as.formula(formula),data=df)
	    hr<-sprintf("%0.3f",exp(coefficients(a)[1]))
	    lci<-sprintf("%0.3f",exp(confint(a)[1,1]))
	    uci<-sprintf("%0.3f",exp(confint(a)[1,2]))
	    hr_ci<- paste0(hr," ","(",lci,"-",uci,")")
	    result<- c(trait,"PRS-CS",phi," ",hr_ci,count,hr)
	    return(result)
	}
	prslist<- c("phi_auto","phi_1","phi_2","phi_3","phi_4")
	training.result<- foreach(i=1:5,.combine="rbind",.packages="survival") %dopar% {
		training_temp<- cox_result(ep,prslist[i],gwas,count,full.training)
	}

	testing.result<- foreach(i=1:5,.combine="rbind",.packages="survival") %dopar% {
		testing_temp<- cox_result(ep,prslist[i],gwas,count,full.testing)
	}
}
