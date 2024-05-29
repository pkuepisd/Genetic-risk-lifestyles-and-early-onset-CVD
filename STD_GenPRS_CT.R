library(data.table)
Args<-commandArgs()
basefile<-Args[3]
clumpfile<-Args[4]
dosagefile<-Args[5]
outfile<-Args[6]

base<-fread(basefile,header=TRUE)
base$EFF<- as.numeric(base$EFF)
base$PVAL<- as.numeric(base$PVAL)

clump<-read.table(clumpfile,as.is=TRUE)
dosage<-fread(dosagefile,header=TRUE)

Dos<- as.matrix(dosage[,7:ncol(dosage)])

colnames(clump)<-"SNP"

baseclump<-subset(base,base$SNP %in% clump$SNP,select=c(SNP,A1,EFF,PVAL))

reverse<-ifelse(baseclump$A1 != dosage$alleleB,TRUE,FALSE)
for (i in 1:nrow(Dos)){
	if (reverse[i]){
		Dos[i,]<-2-Dos[i,]
	}
}

plist <- c(5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 
             seq(0.0001, 0.0009, 0.0001),
             seq(0.001, 0.009, 0.001),
             seq(0.01, 0.09, 0.01),
             seq(0.1, 0.5, 0.1), 1.01)

PRS<-data.frame(id=colnames(dosage)[7:ncol(dosage)])
for (k in plist){
	threshold <- ifelse(baseclump$PVAL<k,TRUE,FALSE)
	effect<-t(as.matrix(baseclump[threshold,"EFF"]))
	df <- Dos[threshold,]
	prs<-as.data.frame(t(effect %*% df))
	PRS<-cbind(PRS,prs)
}
colnames(PRS)[2:ncol(PRS)] <-paste0("prs",seq(1,length(plist)))
fwrite(PRS, file = outfile, row.names = F, sep = " ", quote = F, na = "NA")

