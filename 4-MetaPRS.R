# MetaPRS pipeline
library(data.table)
library(glmnet)
library(doParallel)
library(foreach)
library(survival)
library(plyr)
library(openxlsx)
library(MASS)

stopCluster(cl)
cl<-makeCluster(4)
registerDoParallel(cl)

MetaPRS_calculate<-function(ep,gwaslist,tail,save="FALSE") {
    load(ep,"_full.RData")
    load("training_fam")
    full<- rename(full,c("LDL-C_ALL"="LDL.C_ALL","LDL-C_EAS"="LDL.C_EAS","HDL-C_ALL"="HDL.C_ALL","HDL-C_EAS"="HDL.C_EAS","LDL-C_ALL"="LDL.C_ALL","HDL-C_ALL"="HDL.C_ALL"))
    message("Full data is merged")
    full_training<- subset(full,id %in% training_fam$id)
    full_testing<- subset(full,!id %in% training_fam$id)
    
    sd_trait_ancestry<- NULL
    beta_trait_ancestry<- NULL
    for (k in 1:length(gwaslist)) {
        combine<- strsplit(gwaslist[k],"_")[[1]][2]
        subprs<- strsplit(gwaslist[k],"_")[[1]][1]
        if (combine=="Combine"){
            prs1<- paste0(subprs,"_ALL")
            prs2<- paste0(subprs,"_EAS")
            prsCombine<- paste0(subprs,"_Combine")
            combinelist<- c(prs1,prs2)
            combine_training<- as.matrix(full_training[,c(combinelist,ep)])
            combine_training[,combinelist]<- scale(combine_training[,combinelist])
            sub_sd<- c(sd(full_training[,prs1]),sd(full_training[,prs2]))
            names(sub_sd)<- combinelist
            sd_trait_ancestry<- c(sd_trait_ancestry,sub_sd)
            formula<- paste0(ep,"~",prs1,"+",prs2)
            fullmodel<- glm(formula,family = "binomial",data=as.data.frame(combine_training))
            result<- stepAIC(fullmodel,trace=1)
            final_model<- glm(result$formula,family = "binomial",data=as.data.frame(combine_training))
            sub_beta<- coefficients(final_model)[combinelist]
            sub_beta[is.na(sub_beta)]<- 0
            names(sub_beta)<- combinelist
            beta_trait_ancestry<- c(beta_trait_ancestry,sub_beta)
            newdata<- as.data.frame(sapply(full[,combinelist],scale))
            combinePRS<- predict(final_model,newdata=newdata)
            full[,prsCombine]<- combinePRS
        }
    }
    full_training<- subset(full,id %in% training_fam$id)
    full_testing<- subset(full,!id %in% training_fam$id)
    
    x_training<-as.matrix(full_training[,gwaslist])
    sd_trait<- apply(x_training[,gwaslist],2,sd)
    x_training<- scale(x_training)
    x_testing<-as.matrix(full_testing[,gwaslist])
    x_testing<- scale(x_testing)
    if ("HDL.C_Combine" %in% gwaslist) {
        x_training[,"HDL.C_Combine"]<- -x_training[,"HDL.C_Combine"]
        x_testing[,"HDL.C_Combine"]<- -x_testing[,"HDL.C_Combine"]
    }
    
    y_training<-full_training[,ep]
    type_measure<- "auc"
    
    glmnet.result<- function(alpha){
        set.seed(1) 
        model<-cv.glmnet(x_training,y_training,family="binomial",type.measure = type_measure,nfold=10,parallel = T,alpha=alpha) 
        result<-c(alpha,model[["lambda.min"]],max(model[["cvm"]]))
        return(result)
    }
    
    result<-foreach(i = seq(0,1,0.1),.combine = "rbind", .packages = "glmnet" ) %dopar% {
        result_combine<-glmnet.result(i)
    }
    message("EN training is finished")
    result<-as.data.frame(result)
    colnames(result)<-c("alpha","lambda","auc")
    process<- result
    colnames(process)<- c("Alpha","Lambda","AUC")
    best_alpha<- result[which(result$auc==max(result$auc)),"alpha"]
    best_lambda<- result[which(result$auc==max(result$auc)),"lambda"]
    final_model<-glmnet(x_training,y_training,family="binomial",alpha= best_alpha,lambda = best_lambda)
    coefficients<- as.data.frame(as.matrix(coef(final_model,s=result[which(result$auc==max(result$auc)),"lambda"])))
    beta_trait<- unlist(coefficients)[-1]
    names(beta_trait)<- gwaslist

    coefficients$OR<- exp(coefficients$s1)
   
    MetaPRS_training<- predict(final_model,type="link",newx=x_training)
    MetaPRS_testing<- predict(final_model,type="link",newx=x_testing)
    
    message("MetaPRS is generated")
    
    full_training_meta<- cbind(full_training,MetaPRS_training) 
    full_testing_meta<- cbind(full_testing,MetaPRS_testing) 
    
    colnames(full_training_meta)[length(full_training_meta)]<-"MetaPRS"
    colnames(full_testing_meta)[length(full_testing_meta)]<-"MetaPRS"
    
    prslist<-c("MetaPRS",gwaslist)
    full_training_meta[,prslist]<- scale(full_training_meta[,prslist])
    if ("HDL.C_Combine" %in% prslist) full_training_meta$HDL.C_Combine<- -full_training_meta$HDL.C_Combine
    logistic_result<-function(ep,ex,data){
        fomular<- paste0(ep,"~",ex)
        or_ci<- tryCatch(
            {
                a<- glm(fomular,data=data,family = "binomial")
                or<-sprintf("%0.3f",exp(coefficients(a)[2]))
                lci<-sprintf("%0.3f",exp(confint(a)[2,1]))
                uci<-sprintf("%0.3f",exp(confint(a)[2,2]))
                or_ci<- paste0(or," ","(",lci,"-",uci,")")
            },
            error=function(cond) {
                message(paste0(cond,"\n"))
                return("NA")
            }
        )
        return(or_ci)
    }
    
    result_all<- foreach(i = 1:length(prslist),.combine = "rbind") %dopar% {
        result_combine<-logistic_result(ep,prslist[i],full_training_meta)
    }

    result_all<- as.data.frame(result_all)
    OUT<-cbind(prslist,result_all)
    EN<-sprintf("%0.3f",coefficients$OR)
    EN[1]<-"-"
    OUT<-cbind(OUT,EN)
    colnames(OUT)<-c("PRS","Effect","Elastic Net OR")

    full_meta<- rbind(full_training_meta,full_testing_meta)
    id_ep_meta<- subset(full_meta,select=c("id","MetaPRS",gwaslist))

    if (save==TRUE){
        fwrite(id_ep_meta,paste0("id_",ep,"_ALLMetaPRS_all_",tail,"_binomial.txt"),sep=" ",col.names=TRUE,quote=FALSE)
    }
    weight_info<- list(beta_trait=beta_trait,sd_trait=sd_trait,beta_trait_ancestry=beta_trait_ancestry,sd_trait_ancestry=sd_trait_ancestry)
    output<- list(OUT,weight_info)
    return(output)
}

u1_IHD_amend4<- c("IHD_Combine","IS_Combine","SBP_Combine","DBP_Combine","TC_Combine","TG_Combine","LDL.C_Combine","HDL.C_Combine","FPG_Combine","T2D_Combine")
cad_amend<- MetaPRS_calculate("cad",u1_IHD_amend4,"AM4",save=T)

u1_AIS_amend4<- c("IHD_Combine","IS_Combine","AF_Combine","SBP_Combine","DBP_Combine","TC_Combine","TG_Combine","LDL.C_Combine","HDL.C_Combine","FPG_Combine","T2D_Combine")
is_amend<- MetaPRS_calculate("is",u1_AIS_amend4,"AM4",save=T)

u1_HS_amend4<- c("IHD_Combine","IS_Combine","HS_Combine","SBP_Combine","DBP_Combine","FPG_Combine","T2D_Combine")
hs_amend<- MetaPRS_calculate("hs",u1_HS_amend4,"AM4",save=T)
