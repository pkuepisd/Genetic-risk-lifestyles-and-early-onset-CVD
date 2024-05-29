library(data.table)
Args<-commandArgs()

prsdir<-Args[3]

chr1<-as.data.frame(fread(paste0(prsdir,"/chr1.prs"),header=T))
PRSs<- as.data.frame(chr1[,1])
for (k in 2:41){
	prs<-chr1[,c(1,k)]
	for (i in 2:22){
		df<- fread(paste0(prsdir,"/chr",i,".prs"),header=T)
		df<-as.data.frame(df)
		prs[,i+1]<-df[,k]
	}
	PRSs[,k]<- rowSums(prs[,2:23])
}
colnames(PRSs)<- c("id",paste0("prs",1:40))
write.csv(PRSs, file=paste0(prsdir,"/all.csv"),row.names=F, quote=F)