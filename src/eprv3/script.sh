ulimit -s unlimited
for dfile in {1..6}
do
    for i in {1..10}
    do
        python ppc_eprv3.py -n 1 -dfile $dfile --clust
    done
done