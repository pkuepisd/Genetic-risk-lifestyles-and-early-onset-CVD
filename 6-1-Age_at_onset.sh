# Associations of genetic risk and lifestyles with different age-onset of CVD (sh file) 
eplist="cad is hs"
codepath=~/PhD/code
for ep in $eplist; do
    for index in {1..7};do
    	OUTDIR=~/PhD/Project3/Result/OE
        qsub  $codepath/6-2-Age_at_onset.pbs \
            -N AAO_${ep}_${index} \
            -v ep=$ep,index=$index \
            -l nodes=1:ppn=20 \
            -l walltime=10:00:00 \
            -o $OUTDIR \
            -e $OUTDIR	
    done
done