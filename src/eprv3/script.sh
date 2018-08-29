ulimit -s unlimited
for dfile in {1..3}
do
    for n in {0..3}
    do
        python ppc_eprv3.py -n $n -dfile $dfile --clust
    done
done