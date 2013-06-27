#!/bin/bash
#$ -S /bin/bash
#$ -q bhost.q
source /home/rxzhu/.bashrc
source /home/rxzhu/.bashrc_xiaohan

echo $SHELL

plane=(complex residue)
atg_or_atb=(atg atb)
res_or_tri=(res tri)


for p in "${plane[@]}"
do
    for a in "${atg_or_atb[@]}"
    do
        for r in "${res_or_tri[@]}"
        do
	    echo $p $a $r;
            time /software/python.2.7.3/bin/python $VE_CODE/simmat/fp_175_padded.py $p $a $r > $VE_CODE/simmat/data/"${p}_${a}_${r}.txt";
        done
    done
done
