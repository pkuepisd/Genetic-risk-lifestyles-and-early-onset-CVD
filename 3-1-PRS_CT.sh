# C+T workflow
gwaslist="IHD_ALL IS_ALL HS_ALL AF_ALL SBP_ALL DBP_ALL TC_ALL TG_ALL LDL-C_ALL HDL-C_ALL T2D_ALL FPG_ALL IHD_EAS IS_EAS HS_EAS AF_EAS SBP_EAS DBP_EAS TC_EAS TG_EAS LDL-C_EAS HDL-C_EAS T2D_EAS FPG_EAS"
refdir=~/PhD/Project1/CKB_REF
codepath=~/PhD/code
eplist="cad is hs"
#1. Clumping####################################################################
for gwas in $gwaslist;do
	for qctype in 1;do
		OUTDIR_pre=~/PhD/Project1/GWASdata/$gwas/qctype$qctype
		mkdir $OUTDIR_pre
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2
			mkdir $OUTDIR
			for chr in `seq 22`;do
				gwasfile=~/PhD/Project1/GWASdata/$gwas/$gwas.org.stdqc${qctype}_chr
			    qsub  $codepath/STD_clumping.pbs \
			        -N CT_$gwas.$qctype.$r2.$chr \
			        -v ref=${refdir}/chr${chr},BaseFile=${gwasfile}/chr${chr}.txt,OUT=$OUTDIR/chr${chr},r2=$r2\
			        -l nodes=1:ppn=20 \
			        -l walltime=2:00:00 \
			        -o $OUTDIR \
			        -e $OUTDIR
			done
		done
	done	
done
 
for gwas in $gwaslist;do
	for qctype in 1;do
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2
			for chr in `seq 22`;do
				awk '{print $3}' $OUTDIR/chr${chr}.clumped | tail -n +2 | sed '/^$/d' > $OUTDIR/chr${chr}.snplist
			done
		done
	done	
done

# Calculate PRS###############################################################
for gwas in $gwaslist;do
	BaseDir=~/PhD/Project1/GWASdata/${gwas}/${gwas}.org.stdqc1_chr
	for qctype in 1;do
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2
			for chr in `seq 22`;do
				qsub $codepath/STD_GenPRS_CT.pbs \
			    -N P_$gwas.$qctype.$r2.$chr \
			    -v basefile=$BaseDir/chr${chr}.txt,clumpfile=$OUTDIR/chr${chr}.snplist,dosagefile=$OUTDIR/chr${chr}.dosage,outfile=$OUTDIR/chr${chr}.prs \
			    -l nodes=1:ppn=20 \
			    -l walltime=2:00:00 \
			    -o $OUTDIR \
			    -e $OUTDIR
			done
		done
	done	
done

for gwas in $gwaslist;do
	for qctype in 1;do
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2
			qsub $codepath/STD_SumPRS_CT.pbs \
		    -N S_$gwas.$qctype.$r2.$chr \
		    -v prsdir=$OUTDIR \
		    -l nodes=1:ppn=10 \
		    -l walltime=2:00:00 \
		    -o $OUTDIR \
		    -e $OUTDIR
		done
	done	
done

# Select PRS###########################################################
for gwas in $gwaslist;do
	for qctype in 1;do
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2
			qsub $codepath/Training_PRS_CT_plus.pbs \
		    -N T_$gwas.$qctype.$r2 \
		    -v PRSpath=$OUTDIR,trait=$gwas,r2=$r2,tail=training \
		    -l nodes=1:ppn=20 \
		    -l walltime=2:00:00 \
		    -o $OUTDIR \
		    -e $OUTDIR
		done
	done
done

# Testing set
for gwas in $gwaslist;do
	for r2 in 0 0.2 0.4 0.6 0.8;do
		workdir=~/PhD/Project1/GWASdata/$gwas/qctype1/$r2/testing
		for chr in `seq 22`;do
		    qsub  $codepath/STD_GenPRS_CT_flexible.pbs \
		        -N flexible_$chr \
		        -v gwas=$gwas,chr=$chr,r2=$r2 \
		        -l nodes=1:ppn=20 \
		        -l walltime=20:00:00 \
		        -o $workdir \
		        -e $workdir
		done
	done
done

for gwas in $gwaslist;do
	for qctype in 1;do
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2/testing
			qsub $codepath/STD_SumPRS_CT.pbs \
		    -N S_$gwas.$qctype.$r2.$chr \
		    -v prsdir=$OUTDIR \
		    -l nodes=1:ppn=10 \
		    -l walltime=2:00:00 \
		    -o $OUTDIR \
		    -e $OUTDIR
		done
	done	
done

for gwas in $gwaslist;do
	for qctype in 1;do
		for r2 in 0 0.2 0.4 0.6 0.8;do
			OUTDIR=~/PhD/Project1/GWASdata/$gwas/qctype$qctype/$r2/testing
			qsub $codepath/Training_PRS_CT_plus.pbs \
		    -N T_$gwas.$qctype.$r2 \
		    -v PRSpath=$OUTDIR,trait=$gwas,r2=$r2,tail=testing \
		    -l nodes=1:ppn=20 \
		    -l walltime=2:00:00 \
		    -o $OUTDIR \
		    -e $OUTDIR
		done
	done	
done