for i in {1..10}
do
    python -m cProfile -o profile.txt ppc_eprv3.py -n 2 -nlive 25 -nrep 3 -prec 0.01
done

python profiler.py