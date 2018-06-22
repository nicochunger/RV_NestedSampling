for i in {1..2}
do
    #python ppc_eprv3.py
    python -m cProfile -o profile.txt ppc_eprv3.py
done

python profiler.py