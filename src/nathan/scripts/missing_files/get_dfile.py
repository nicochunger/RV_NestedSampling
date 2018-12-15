import sys
import subprocess

data_path = '/home/nunger/tesis/codigo/src/nathan/data/data_2gen'
# data_path = '/home/spectro/nunger/codigo/src/nathan/data/data_2gen'
data_files = subprocess.check_output(
    'ls {}'.format(data_path), shell=True).decode('utf-8').split('\n')
data_files.remove('')
data_files.remove('prueba.txt')

res_path = '/home/nunger/tesis/resultados/nathan/results_2gen'
# res_path = '/scratch/nunger/nathan/results_2gen'
res_files = subprocess.check_output(
    'ls {}'.format(res_path), shell=True).decode('utf-8').split('\n')
res_files.remove('')
# res_files.remove('results_2gen.zip')

real_res = []
for dfile in data_files:
    for model in range(1, 5):
        real_res.append('results_'+dfile[5:-4]+'_model{}'.format(model)+'.txt')

# Find missing result files
missing = []
for d in real_res:
    if d not in res_files:
        missing.append(d)

idxs = []
for mis in missing:
    dfile = 'data_'+mis[8:-11]+'.txt'
    idxs.append((data_files.index(dfile)+1, int(mis[-5])))

print(idxs)

# # Get datafile and option
# # opt = 0 is the dfile
# # opt = 1 is the model
# idx = int(sys.argv[1])
# opt = int(sys.argv[2])

# # Print requested result for bash output
# print(idxs[idx][opt])
