# PGS catalog pipeline
codepath=~/PhD/code
testing=~/PhD/Project1/testing.fam
PGSlist="PGS003725 PGS000337 PGS002262 PGS000013 PGS002259 PGS000039"

for pgs in $PGSlist; do
    qsub ~/PhD/Project1/code/PGSFileStandardQC.pbs \
        -v pgs=$pgs,INFO=0.8,FREQ=0.01,HWP=1e-6 \
        -N $pgs \
        -l nodes=1:ppn=10 \
        -l walltime=1:00:00
done

for PGS in $PGSlist;do
	OUTDIR=~/PhD/Project1/PGSdata/$PGS.txt.org.stdqc_chr
	cd $OUTDIR
	chrs=`ls *.txt | sed "s/.txt//g" | sed "s/chr//g" | sort -n`
	cd
	for chr in $chrs;do
		qsub $codepath/STD_GenPRS_DirLD_flexible.pbs \
	    -N P_$PGS.$chr \
	    -v PGS=$PGS,chr=$chr \
	    -l nodes=1:ppn=10 \
	    -l walltime=3:00:00 \
	    -o $OUTDIR \
	    -e $OUTDIR
	done
done

for PGS in $PGSlist;do
	OUTDIR=~/PhD/Project1/PGSdata/$PGS.txt.org.stdqc_chr
	qsub $codepath/STD_SumPRS_DirLD.pbs \
    -N S_$PGS \
    -v PGS=$PGS \
    -l nodes=1:ppn=5 \
    -l walltime=3:00:00 \
    -o $OUTDIR \
    -e $OUTDIR
done

R --no-save <$codepath/Training_PRS_DirLD.R