library(data.table)
library(plyr)
Args<-commandArgs()

PGS<- Args[3]

prsdir<-paste0("~/PhD/Project1/PGSdata/",PGS,".txt.org.stdqc_chr/")
chrname<- list.files(path = prsdir, pattern = "*.prs")

for (i in 1:length(chrname)){
	df<- fread(paste0(prsdir,chrname[i]),header=T)
	colnames(df)[2]<- chrname[i]
	if (i==1) {
		PRS<- df
	} else {
		PRS<- merge(PRS,df,by="id")
	}
}
PRS$V1 <- rowSums(PRS[,2:ncol(PRS)])
PRS<- subset(PRS,select=c(id,V1))
colnames(PRS)<- c("id",PGS)
write.csv(PRS, file=paste0(prsdir,"/all.csv"), row.names=F, quote=F)

