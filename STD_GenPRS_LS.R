library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)
library(doParallel)
library(survival)
cl <- makeCluster(10)
registerDoParallel(cl)

Args<- commandArgs()
gwas<- Args[3]

basedir <- paste0("~/PhD/Project1/GWASdata/",gwas,"/",gwas,".org.stdqc1_chr")
testdir<- "/public/data/GWAS/2018-11-22-Impute_QCed/gwasplink"

if (strsplit(gwas,"_")[[1]][2]=="ALL" | strsplit(gwas,"_")[[1]][2]=="ALL2"){
	LDblocks <- "EUR.hg19"
	refdir<- "/public/data/1000G_EUR"
} else {
	LDblocks <- "ASN.hg19"
	refdir<- "/public/home/yangsongchun/1000genomes/phase3v5a/EAS"
}

gwas_samplesize<- read.csv("~/PhD/Project1/GWAS_samplesize.csv",header=T)
n<- gwas_samplesize[gwas_samplesize$gwaslist==gwas,"samplesize"]

# Gen PRS
for (i in 1:22) {
	message(paste0("Chr ",i," is calculating"))
	ss<-fread(paste0(basedir,"/chr",i,".txt"))
	ss$PVAL<- as.numeric(ss$PVAL)
	ref.bfile<- paste0(refdir,"/chr",i)
	ref.bim<- fread(paste0(refdir,"/chr",i,".bim"))
	test.bfile<- paste0(testdir,"/chr",i)
	if (length(table(ss$PVAL==0))>1) {
		ss$tag<- 1:nrow(ss)
		ss_p0<- ss[ss$PVAL==0,]
		ss_p1<- ss[ss$PVAL!=0,]
		ss_p1$cor<- p2cor(ss_p1$PVAL,n,ss_p1$EFF)
		t_p0<- ss_p0$EFF/ss_p0$SE
		ss_p0$cor<- t_p0/sqrt(n-2+t_p0^2)
		ss<- rbind(ss_p0,ss_p1)
		ss<- ss[order(ss$tag),]
		cor<- ss$cor
	} else {
		cor <- p2cor(p = ss$PVAL, n = n, sign=ss$EFF)	
	}

	if (i==1) {
		out <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$POS, 
                     A1=ss$A1, A2=ss$A2, 
                     ref.bfile=ref.bfile, test.bfile=test.bfile,
                     LDblocks = LDblocks,cluster=cl)
	} else {
		out_temp <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$POS, 
                     A1=ss$A1, A2=ss$A2,
                     ref.bfile=ref.bfile, test.bfile=test.bfile,
                     LDblocks = LDblocks,cluster=cl)
		message("Prepare to merge")
		out<- merge(out,out_temp)
	}
}