# QC of the GWAS data (sh file)
gwaslist="IHD_ALL IS_ALL HS_ALL AF_ALL SBP_ALL DBP_ALL TC_ALL TG_ALL LDL-C_ALL HDL-C_ALL T2D_ALL FPG_ALL IHD_EAS IS_EAS HS_EAS AF_EAS SBP_EAS DBP_EAS TC_EAS TG_EAS LDL-C_EAS HDL-C_EAS T2D_EAS FPG_EAS"
target1=~/PhD/Project1/ckb_info_table_modified_V3.txt
codepath=~/PhD/code 

for gwas in $gwaslist; do
    gwasfile=~/PhD/Project1/GWASdata/$gwas/$gwas.org
    qsub $codepath/2-2-STD_GWASQC.pbs \
        -v gwas=$gwas,gwasfile=$gwasfile,ssfinfo=0.8,ssfmaf=0.01,target=$target1,da="TRUE",di="TRUE",ckinfo="TRUE",ckinfo2=0.8,ckfreq="TRUE",ckfreq2=0.01,ckhwe="TRUE",ckhwe2=1e-6,qctype=1 \
        -N $gwas.1 \
        -l nodes=1:ppn=6 \
        -l walltime=1:00:00 \
        -o ~/PhD/Project1/GWASdata/${gwas} \
        -e ~/PhD/Project1/GWASdata/${gwas}
done