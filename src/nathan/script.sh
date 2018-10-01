ulimit -s unlimited
for dfile in {4..6}
do
    for n in {1..6}
    do
        echo Datafile: $dfile, Repetition: $n/6
        python ppc_nathan.py -model 2 -dfile $dfile
    done
done