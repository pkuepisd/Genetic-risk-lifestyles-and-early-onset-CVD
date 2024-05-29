#!bin/bash
for ep in cad is hs;do
    qsub  /public/home/sundong/PhD/code/8-2-IR_adj_boot.pbs \
        -N G10-1_${ep} \
        -v ep=$ep \
        -l nodes=1:ppn=20 \
        -l walltime=10:00:00 \
        -o $OUTDIR \
        -e $OUTDIR
done