ulimit -s unlimited
for dfile in {2..2}
do
    for n in {0..5}
    do
        python ppc_eprv3.py -n 1 -dfile $dfile --clust
    done
done