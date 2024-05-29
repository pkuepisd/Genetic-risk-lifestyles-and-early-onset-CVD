# Compare trait-specific PRS
library(data.table)
library(doParallel)
library(foreach)
library(survival)
library(ggplot2)

stopCluster(cl)
cl<- makeCluster(4)
registerDoParallel(cl)

cox_result<- function(ep,ex,data){
    formula<- paste0("Surv(age_enter,age_",ep,",",ep,")~",ex,"+national_pc01+national_pc02+national_pc03+national_pc04+national_pc05+national_pc06+national_pc07+national_pc08+national_pc09+national_pc10+array+education_3g+marital_2g+is_female")
    a<-coxph(as.formula(formula),data=data)
    hr<- exp(coefficients(a)[1])
    lci<- exp(confint(a)[1,1])
    uci<- exp(confint(a)[1,2])
    hr_ci<- paste0(sprintf("%0.2f",hr)," ","(",sprintf("%0.2f",lci),"-",sprintf("%0.2f",uci),")")
    result<- c(ep,ex,hr_ci,hr,lci,uci)
    return(result)
}

for (ep in c("cad","is","hs")){
    load(paste0(ep,"_prsfile.RData"))
    load(paste0(ep,"_full_testing.RData"))
    gwaslist<- colnames(prsfile)[-1]
    
    result<- foreach (i=1:length(gwaslist),.combine="rbind",.packages="survival") %dopar% {
        bind<- cox_result(ep,gwaslist[i],full_testing)
    }
}

