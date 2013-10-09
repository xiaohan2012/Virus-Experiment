#types=(spe spa spb tma tmb tmc)
dirname=$1

shift

for t in "$@"
do
    /software/python.2.7.3/bin/python $VE_CODE/simmat/source.py "data/$dirname/$t.txt" > "data/$dirname/$t.csv"
done
