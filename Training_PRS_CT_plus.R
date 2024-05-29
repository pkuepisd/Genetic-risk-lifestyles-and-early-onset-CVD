library(data.table)
library(doParallel)
library(foreach)
library(survival)
cl<-makeCluster(20)
registerDoParallel(cl)
cox_result<- function(ep,ex,trait,r2,p){
	formula<- paste0("Surv(age_enter,age_",ep,",",ep,")~",ex,"+national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array+strata(is_female)")
    a<-coxph(as.formula(formula),data=training_PRS)
    hr<-sprintf("%0.3f",exp(coefficients(a)[1]))
    lci<-sprintf("%0.3f",exp(confint(a)[1,1]))
    uci<-sprintf("%0.3f",exp(confint(a)[1,2]))
    hr_ci<- paste0(hr," ","(",lci,"-",uci,")")
    result<- c(trait,"C+T",r2,p,hr_ci)
    return(result)
}
Args<- commandArgs()

PRSpath<- Args[3]
trait<- Args[4]
r2<- Args[5]
tail<- Args[6]

load("training_PRS.RData")

eplist<-c("cad","is","hs")
plist<- c(5e-08,1e-07,5e-07,1e-06,5e-06,1e-05,5e-05,1e-04,2e-04,3e-04,4e-04,5e-04,6e-04,7e-04,8e-04,9e-04,1e-03,2e-03,3e-03,4e-03,5e-03,6e-03,7e-03,8e-03,9e-03,1e-02,2e-02,3e-02,4e-02,5e-02,6e-02,7e-02,8e-02,9e-02,1e-01,2e-01,3e-01,4e-01,5e-01,1)

for (EP in eplist) {
	result<- foreach(i=1:40,.combine="rbind",.packages="survival") %dopar% {
		result_combine<- cox_result(EP,paste0("prs",i),trait,r2,plist[i])
	}
    colnames(result)<- c("Trait","Method","Parameter1","Parameter2","HR_CI")
	fwrite(result,paste0(PRSpath,"/CT_",EP,"_",tail,".csv"),row.names = F, sep = ",", quote = F, na = "NA")
}
stopCluster(cl)