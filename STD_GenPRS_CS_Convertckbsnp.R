library(data.table)
Args<- commandArgs()

gwas<- Args[3]

CSpath<- paste0("~/PhD/Project1/GWASdata/",gwas,"/",gwas,".org.stdqc1_chr/CS/")
info<- fread("~/PhD/Project1/ckb_info_table_modified_V3.txt",header=TRUE)

for (chr in 1:22) {
	message(paste0("Chr",chr," is converting"))
	for (phi in c("auto","1","1e-2","1e-4","1e-6")){
		message(paste0("Phi",phi," is converting"))

		if (phi=="auto") {
			CS<- fread(paste0(CSpath,"chr",chr,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"),header=FALSE)
		}
		if (phi=="1") {
			CS<- fread(paste0(CSpath,"chr",chr,"_1_pst_eff_a1_b0.5_phi1e+00_chr",chr,".txt"),header=FALSE)
		}
		if (phi=="1e-2") {
			CS<- fread(paste0(CSpath,"chr",chr,"_1e-2_pst_eff_a1_b0.5_phi1e-02_chr",chr,".txt"),header=FALSE)
		}
		if (phi=="1e-4") {
			CS<- fread(paste0(CSpath,"chr",chr,"_1e-4_pst_eff_a1_b0.5_phi1e-04_chr",chr,".txt"),header=FALSE)
		}
		if (phi=="1e-6") {
			CS<- fread(paste0(CSpath,"chr",chr,"_1e-6_pst_eff_a1_b0.5_phi1e-06_chr",chr,".txt"),header=FALSE)
		}

		colnames(CS)<- c("CHR","rsid","BP","A1","A2","A1_Effect")


		if (phi=="auto"){
			CS$tag<- 1:nrow(CS)
			sub_info<- subset(info,CHR=chr,select=c("rsid","SNP"))
			CS_info<- merge(CS,sub_info,by="rsid",all.x=TRUE)
			CS_info<- CS_info[order(CS_info$tag)]
			SNP<- CS_info$SNP
			
		} else {
			CS_info<- cbind(CS,SNP)
		}

		final<- CS_info[,c("CHR","SNP","BP","A1","A2","A1_Effect")]
		fwrite(final,paste0(CSpath,"chr",chr,"_",phi,".txt"),col.names=FALSE,sep=" ",quote=FALSE)
	}
}