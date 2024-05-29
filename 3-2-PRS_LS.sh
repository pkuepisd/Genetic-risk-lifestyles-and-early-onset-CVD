# lassosum workflow
gwaslist="IHD_ALL IS_ALL HS_ALL AF_ALL SBP_ALL DBP_ALL TC_ALL TG_ALL LDL-C_ALL HDL-C_ALL T2D_ALL FPG_ALL IHD_EAS IS_EAS HS_EAS AF_EAS SBP_EAS DBP_EAS TC_EAS TG_EAS LDL-C_EAS HDL-C_EAS T2D_EAS FPG_EAS"
eplist="cad is hs"

for gwas in $gwaslist;do
	OUTDIR=~/PhD/Project1/GWASdata/$gwas/$gwas.org.stdqc1_chr
    qsub  $codepath/STD_GenPRS_LS.pbs \
        -N LS_$gwas \
        -v gwas=$gwas \
        -l nodes=1:ppn=20 \
        -l walltime=100:00:00 \
        -o $OUTDIR \
        -e $OUTDIR	
done

for gwas in $gwaslist;do
    for ep in $eplist;do
        OUTDIR=~/PhD/Project1/GWASdata/$gwas/$gwas.org.stdqc1_chr
        qsub  $codepath/Training_PRS_LS.pbs \
            -N T_LS_${gwas}_$ep \
            -v gwas=$gwas,ep=$ep \
            -l nodes=1:ppn=20 \
            -l walltime=10:00:00 \
            -o $OUTDIR \
            -e $OUTDIR  
    done
done