library(data.table)
library(survival)

cox_result<- function(ep,ex,method){
	formula<- paste0("Surv(age_enter,age_",ep,",",ep,")~",ex,"+national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array+is_female+education_3g+marital_2g")
    a<-coxph(as.formula(formula),data=testing_PRS)
    hr<-sprintf("%0.3f",exp(coefficients(a)[1]))
    lci<-sprintf("%0.3f",exp(confint(a)[1,1]))
    uci<-sprintf("%0.3f",exp(confint(a)[1,2]))
    hr_ci<- paste0(hr," ","(",lci,"-",uci,")")
    result<- c(ex,method," "," ",hr_ci,hr)
    return(result)
}

load("testing_PRS.RData")
eplist<-c("cad","is","hs")
PGSlist<- c("PGS003725", "PGS000337", "PGS002262", "PGS000013", "PGS002259", "PGS000039")

for (EP in eplist) {
    result1<- cox_result(EP,PGSlist[1],"MetaPRS")
    result2<- cox_result(EP,PGSlist[2],"C+T")
    result3<- cox_result(EP,PGSlist[3],"MetaPRS")
    result4<- cox_result(EP,PGSlist[4],"LDpred")
    result5<- cox_result(EP,PGSlist[5],"MetaPRS")
    result6<- cox_result(EP,PGSlist[6],"MetaPRS")
}