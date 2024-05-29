# PRS-CS pipeline
gwaslist="IHD_ALL IS_ALL HS_ALL AF_ALL SBP_ALL DBP_ALL TC_ALL TG_ALL LDL-C_ALL HDL-C_ALL T2D_ALL FPG_ALL IHD_EAS IS_EAS HS_EAS AF_EAS SBP_EAS DBP_EAS TC_EAS TG_EAS LDL-C_EAS HDL-C_EAS T2D_EAS FPG_EAS"
eplist="cad is hs"
refdir_ALL=~/PhD/Project1/PRSCS/PRScs/ldblk_1kg_eur
refdir_EAS=~/PhD/Project1/PRSCS/PRScs/ldblk_1kg_eas
bimdir=~/PhD/Project1/PRSCS/newbim
pydir=~/PhD/Project1/PRSCS/PRScs
gwaslist_ALL=("IHD_ALL" "IS_ALL" "HS_ALL" "AF_ALL" "SBP_ALL" "DBP_ALL" "TG_ALL" "TC_ALL" "LDL-C_ALL" "HDL-C_ALL" "T2D_ALL" "FPG_ALL")
Nlist_ALL=("547261" "1296908" "456348" "1030836" "757601" "757601" "1320016" "1320016" "1320016" "1320016" "898130" "200622")
gwaslist_EAS=("IHD_EAS" "IS_EAS" "HS_EAS" "AF_EAS" "SBP_EAS" "DBP_EAS" "TG_EAS" "TC_EAS" "LDL-C_EAS" "HDL-C_EAS" "T2D_EAS" "FPG_EAS")
Nlist_EAS=("168228" "174686" "153478" "36792" "145505" "145515" "111667" "135808" "72866" "74970" "210865" "288127")

# _ALL
for i in {0..11};do
	gwas=${gwaslist_ALL[i]}
	N=${Nlist_ALL[i]}
	ssfdir=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
	OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
	for chr in `seq 22`;do
	    qsub  $codepath/STD_GenPRS_CS.pbs \
	        -N CS_${gwas}_${chr} \
	        -v pydir=$pydir,refdir=$refdir_ALL,bimdir=$bimdir,chr=$chr,ssfdir=$ssfdir,N=$N \
	        -l nodes=1:ppn=20 \
	        -l walltime=20:00:00 \
	        -o $OUTDIR \
	        -e $OUTDIR	
    done
done
# _EAS
for i in {0..11};do
	gwas=${gwaslist_EAS[i]}
	N=${Nlist_EAS[i]}
	ssfdir=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
	OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
	for chr in `seq 22`;do
	    qsub  $codepath/STD_GenPRS_CS.pbs \
	        -N CS_${gwas}_${chr} \
	        -v pydir=$pydir,refdir=$refdir_EAS,bimdir=$bimdir,chr=$chr,ssfdir=$ssfdir,N=$N \
	        -l nodes=1:ppn=20 \
	        -l walltime=20:00:00 \
	        -o $OUTDIR \
	        -e $OUTDIR	
    done
done

# grid search:1 1e-2 1e-4 1e-6
for phi in 1 1e-2 1e-4 1e-6;do
	for i in {0..11};do
		gwas=${gwaslist_ALL[i]}
		N=${Nlist_ALL[i]}
		ssfdir=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
		OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
		for chr in `seq 22`;do
		    qsub  $codepath/STD_GenPRS_CS_grid_search.pbs \
		        -N CS_${gwas}_chr${chr}_phi${phi} \
		        -v pydir=$pydir,refdir=$refdir_ALL,bimdir=$bimdir,chr=$chr,ssfdir=$ssfdir,N=$N,phi=$phi \
		        -l nodes=1:ppn=20 \
		        -l walltime=20:00:00 \
		        -o $OUTDIR \
		        -e $OUTDIR	
	    done
	done
done


for phi in 1 1e-2 1e-4 1e-6;do
	for i in {0..11};do
		gwas=${gwaslist_EAS[i]}
		N=${Nlist_EAS[i]}
		ssfdir=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
		OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
		for chr in `seq 22`;do
		    qsub  $codepath/STD_GenPRS_CS_grid_search.pbs \
		        -N CS_${gwas}_chr${chr}_phi${phi} \
		        -v pydir=$pydir,refdir=$refdir_EAS,bimdir=$bimdir,chr=$chr,ssfdir=$ssfdir,N=$N,phi=$phi \
		        -l nodes=1:ppn=20 \
		        -l walltime=20:00:00 \
		        -o $OUTDIR \
		        -e $OUTDIR	
	    done
	done
done

for gwas in $gwaslist;do
	OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
    qsub  $codepath/STD_GenPRS_CS_Convertckbsnp.pbs \
        -N CSconvert_$gwas \
        -v gwas=$gwas \
        -l nodes=1:ppn=10 \
        -l walltime=5:00:00 \
        -o $OUTDIR \
        -e $OUTDIR	
done

for gwas in $gwaslist;do
	for CHR in `seq 22`;do
		for phi in auto 1 1e-2 1e-4 1e-6;do
			OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
		    qsub  $codepath/STD_GenPRS_CS_Genchrprs.pbs \
		        -N CSprs_$gwas.$CHR.$phi \
		        -v gwas=$gwas,CHR=$CHR,phi=$phi \
		        -l nodes=1:ppn=20 \
		        -l walltime=5:00:00 \
		        -o $OUTDIR \
		        -e $OUTDIR
	    done
	done
done

for gwas in $gwaslist;do
	OUTDIR=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr/CS
    qsub  $codepath/STD_SumPRS_CS.pbs \
        -N Sumprs_$gwas \
        -v gwas=$gwas \
        -l nodes=1:ppn=20 \
        -l walltime=5:00:00 \
        -o $OUTDIR \
        -e $OUTDIR
done