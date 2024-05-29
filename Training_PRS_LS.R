library(data.table)
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

Args<- commandArgs()
gwas<- Args[3]
ep<- Args[4]

load(gwas,"_full.training.RData")
load(gwas,"_full.testing.RData")

prslist<- c("X0.2.1","X0.2.2","X0.2.3","X0.2.4","X0.2.5","X0.2.6","X0.2.7","X0.2.8","X0.2.9","X0.2.10","X0.2.11","X0.2.12","X0.2.13","X0.2.14","X0.2.15","X0.2.16","X0.2.17","X0.2.18","X0.2.19","X0.2.20","X0.5.1","X0.5.2","X0.5.3","X0.5.4","X0.5.5","X0.5.6","X0.5.7","X0.5.8","X0.5.9","X0.5.10","X0.5.11","X0.5.12","X0.5.13","X0.5.14","X0.5.15","X0.5.16","X0.5.17","X0.5.18","X0.5.19","X0.5.20","X0.9.1","X0.9.2","X0.9.3","X0.9.4","X0.9.5","X0.9.6","X0.9.7","X0.9.8","X0.9.9","X0.9.10","X0.9.11","X0.9.12","X0.9.13","X0.9.14","X0.9.15","X0.9.16","X0.9.17","X0.9.18","X0.9.19","X0.9.20","X1.1","X1.2","X1.3","X1.4","X1.5","X1.6","X1.7","X1.8","X1.9","X1.10","X1.11","X1.12","X1.13","X1.14","X1.15","X1.16","X1.17","X1.18","X1.19","X1.20")

cox_result<- function(ep,ex,trait,index,df){
	if (index<=20) {
		s<- "0.2"
		col<- index
	}
	if (index>20 & index <=40) {
		s<- "0.5"
		col<- index-20
	}
	if (index>40 & index <=60) {
		s<- "0.9"
		col<- index-40
	}
	if (index>60 & index <=80) {
		s<- "1"
		col<- index-60
	}
	lambda.test<- lambda[col]
	summary<- weight_adj[[paste0(s,".",index)]]
	count<- length(summary[summary!=0])
    formula<- paste0("Surv(age_enter,age_",ep,",",ep,")~",ex,"+national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array+strata(is_female)")
    a<-coxph(as.formula(formula),data=df)
    hr<-sprintf("%0.3f",exp(coefficients(a)[1]))
    lci<-sprintf("%0.3f",exp(confint(a)[1,1]))
    uci<-sprintf("%0.3f",exp(confint(a)[1,2]))
    hr_ci<- paste0(hr," ","(",lci,"-",uci,")")
    result<- c(trait,"lassosum",s,lambda.test,hr_ci,count,hr)
    return(result)
}

training.result<- foreach(i=1:80,.combine="rbind",.packages="survival") %dopar% {
	training_temp<- cox_result(ep,prslist[i],gwas,i,full.training)
}

testing.result<- foreach(i=1:80,.combine="rbind",.packages="survival") %dopar% {
	testing_temp<- cox_result(ep,prslist[i],gwas,i,full.testing)
}
