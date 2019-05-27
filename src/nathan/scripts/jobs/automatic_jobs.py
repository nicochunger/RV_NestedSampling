# Script that automatically determines which datasets and models have to be run

import subprocess
import numpy as np
import sys
import os

# dfile = int(sys.argv[1])
# model = int(sys.argv[2])


# create datafile list and number
# data_path = "/home/nunger/tesis/codigo/src/nathan/data/data_2gen"
data_path = "~/codigo/src/nathan/data/data_2gen"
data_files = subprocess.check_output(
    'ls {}'.format(data_path), shell=True).decode('utf-8').split('\n')
data_files.remove('')  # Remove last empty item
# clean up
for i, dat in enumerate(data_files):
    data_files[i] = data_files[i].strip("data_")
    data_files[i] = data_files[i].strip(".txt")


# data structure to save which datafiles have been run, for what models and how many times
ndfiles = 200
nmodels = 4
jobs_done = np.zeros([ndfiles, nmodels])

# parse results folder
# results_path = "/home/nunger/tesis/resultados/nathan/results_lognormal"
results_path = "/scratch/nunger/nathan/results_jeff-16"
try:
    jobs = subprocess.check_output(
        "wc -l {}/results*".format(results_path), shell=True).decode("utf-8").split("\n")
    # delete last two lines
    jobs.pop(-1)
    jobs.pop(-1)

    for i in range(len(jobs)):
        jobs[i] = jobs[i].strip()

    for i, data in enumerate(data_files):
        for job in jobs:
            if data in job:
                model = int(job[-5])-1
                jobs_done[i, model] = int(job.split(" ")[0])-1
except:
    pass

for i in range(ndfiles):
    print("Datafile "+str(i+1)+": "+str(jobs_done[i, :]))
# print(jobs_done)

# Send missing jobs

for dfile in range(1):
    for model in range(nmodels):
        print("Datafile: "+str(dfile+1)+"; Model: "+str(model+1))
        # Check is jobs are missing
        if jobs_done[dfile, model] < 7:
            # Edit job file for corresponding dfile and model
            with open("ciecc", "r") as f:
                lines = f.readlines()
                splitted = lines[-1].split(" ")
                splitted[5] = str(dfile+1)
                splitted[7] = str(model+1)
                joined = " ".join(splitted)
                lines[-1] = joined

            with open("ciecc", "w") as f:
                f.write("".join(lines))

            nmissing = int(7 - jobs_done[dfile, model])
            # if nmissing < 0:
            #     nmissing = 0
            for _ in range(nmissing):
                os.system("sbatch ciecc")
# sys.exit(nmissing)
