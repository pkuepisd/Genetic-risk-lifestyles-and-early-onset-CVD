library(data.table)
Args<-commandArgs()

gwas<- Args[3]
chr<-Args[4]
r2<- as.numeric(Args[5])

basefile<-paste0("~/PhD/Project1/GWASdata/",gwas,"/",gwas,".org.stdqc1_chr/chr",chr,".txt")
clumpfile<-paste0("~/PhD/Project1/GWASdata/",gwas,"/qctype1/",r2,"/chr",chr,".snplist")
outfile<-paste0("~/PhD/Project1/GWASdata/",gwas,"/qctype1/",r2,"/testing/chr",chr,".prs")

splitpath<-paste0("~/PhD/Project1/GWASdata/",gwas,"/qctype1/",r2,"/testing/chr",chr,".dosage.split")
files <- list.files(path = splitpath, pattern = "sub*")

idfile<- paste0("~/PhD/Project1/GWASdata/",gwas,"/qctype1/",r2,"/testing/chr",chr,".dosage.head")
head <- fread(idfile, header=F)
id <- as.numeric(head[,7:ncol(head)])
PRS<-data.frame(id)

base<-fread(basefile,header=TRUE)
base$EFF<- as.numeric(base$EFF)
base$PVAL<- as.numeric(base$PVAL) 
clump<-read.table(clumpfile,as.is=TRUE)
colnames(clump)<-"SNP"
baseclump_all<-subset(base,base$SNP %in% clump$SNP,select=c(SNP,A1,EFF,PVAL))

i<- 1
for (sub in files) {
	message(paste0("sub= ",sub))
	message(paste0("i= ",i))
	dosagefile<-paste0("~/PhD/Project1/GWASdata/",gwas,"/qctype1/",r2,"/testing/chr",chr,".dosage.split/",sub)
	dosage<-fread(dosagefile,header=FALSE)
	colnames(dosage)[1:6]<-c("chromosome","SNPID","rsid","position","alleleA","alleleB")
	Dos<- as.matrix(dosage[,7:ncol(dosage)])

	baseclump<- baseclump_all[baseclump_all$SNP %in% dosage$rsid,]
	message(nrow(baseclump))
	reverse<-ifelse(baseclump$A1 != dosage$alleleB,TRUE,FALSE)
	for (j in 1:nrow(Dos)){
		if (reverse[j]){
			Dos[j,]<-2-Dos[j,]
		}
	}
	plist <- c(5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 
             seq(0.0001, 0.0009, 0.0001),
             seq(0.001, 0.009, 0.001),
             seq(0.01, 0.09, 0.01),
             seq(0.1, 0.5, 0.1), 1.01)

	for (k in 1:length(plist)) {
		message(paste0("k= ",k))
		threshold <- ifelse(baseclump$PVAL<plist[k],TRUE,FALSE)
		effect<-t(as.matrix(baseclump[threshold,"EFF"]))
		df <- Dos[threshold,]
		prs<-as.data.frame(t(effect %*% df))
		if (i==1) {
			PRS<-cbind(PRS,prs)
		} else {
			PRS[,k+1]<- PRS[,k+1]+prs
		}
	}
	i<- i+1
	rm(dosage)
	rm(Dos)
	rm(df)
	gc()
}
colnames(PRS)[2:ncol(PRS)] <-paste0("prs",seq(1,length(plist)))
fwrite(PRS, file = outfile, row.names = F, sep = " ", quote = F, na = "NA")
