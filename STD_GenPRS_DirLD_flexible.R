library(data.table)
Args<-commandArgs()

PGS<- Args[3]
chr<-Args[4]

basefile<-paste0("~/PhD/Project1/PGSdata/",PGS,".txt.org.stdqc_chr/chr",chr,".txt")
outfile<-paste0("~/PhD/Project1/PGSdata/",PGS,".txt.org.stdqc_chr/chr",chr,".prs")

splitpath<-paste0("~/PhD/Project1/PGSdata/",PGS,".txt.org.stdqc_chr/chr",chr,".dosage.split")
files <- list.files(path = splitpath, pattern = "sub*")

idfile<- paste0("~/PhD/Project1/PGSdata/",PGS,".txt.org.stdqc_chr/chr",chr,".dosage.head")
head <- fread(idfile, header=F)
id <- as.numeric(head[,7:ncol(head)])
PRS<-data.frame(id)

base<-fread(basefile,header=TRUE)
base$BETA<- as.numeric(base$BETA)

i<- 1
for (sub in files) {
	message(paste0("sub= ",sub))
	message(paste0("i= ",i))
	dosagefile<-paste0(splitpath,"/",sub)
	dosage<-fread(dosagefile,header=FALSE)
	colnames(dosage)[1:6]<-c("chromosome","SNPID","rsid","position","alleleA","alleleB")
	Dos<- as.matrix(dosage[,7:ncol(dosage)])
	base_select<- base[base$SNP %in% dosage$rsid,]
	message(nrow(base_select))
	reverse<-ifelse(base_select$EA != dosage$alleleB,TRUE,FALSE)
	for (j in 1:nrow(Dos)){
		if (reverse[j]){
			Dos[j,]<- 2-Dos[j,]
		}
	}
	Dos[is.na(Dos)] <- 0
	effect<-t(as.matrix(base_select$BETA))
	prs<-as.data.frame(t(effect %*% Dos))

	if (i==1) {
		PRS<-cbind(PRS,prs)
	} else {
		PRS[,2]<- PRS[,2]+prs
	}
	i<- i+1
	rm(dosage)
	rm(Dos)
	gc()
}

fwrite(PRS, file = outfile, row.names = F, sep = " ", quote = F, na = "NA")
