plane=(complex residue)
atg_or_atb=(atg atb)
res_or_tri=(res tri)


for p in "${plane[@]}"
do
    for a in "${atg_or_atb[@]}"
    do
        for r in "${res_or_tri[@]}"
        do
	    filename="${p}_${a}_${r}"
            #/software/python.2.7.3/bin/python $VE_CODE/simmat/source.py $VE_CODE/simmat/data/$filename.txt > $VE_CODE/simmat/data/fp_175_padded/$filename.csv;
            /software/python.2.7.3/bin/python $VE_CODE/simmat/fp_175_padded.py print_fp $p $a $r > $VE_CODE/simmat/data/fp_175_padded/$filename.fp.csv;

        done
    done
done
