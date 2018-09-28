ulimit -s unlimited
for dfile in {1..6}
do
    for n in {1..3}
    do
        python ppc_eprv3.py -n 1 -dfile $dfile --narrow
    done
done