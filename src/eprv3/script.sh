ulimit -s unlimited
for dfile in {1..3}
do
    for i in {1..3}
    do
        python ppc_eprv3.py -n 3 -nlive 25 -dfile $dfile --clust
    done
done