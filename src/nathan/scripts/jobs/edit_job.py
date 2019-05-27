# Takes ciecc job file and edits the file for other datafiles/models

import sys

dfile = sys.argv[1]
model = sys.argv[2]

# Open file
with open("ciecc", "r") as f:
    lines = f.readlines()
    splitted = lines[-1].split(" ")
    splitted[5] = dfile
    splitted[7] = model
    joined = " ".join(splitted)
    lines[-1] = joined

with open("ciecc", "w") as f:
    f.write("".join(lines))
