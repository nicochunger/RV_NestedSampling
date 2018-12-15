import subprocess

misfiles = 145
scripts = 15
files_per_script = int(misfiles/scripts)
files_last_scripts = misfiles % scripts

# Write bash scripts
for i in range(scripts+1):
    file = open('mis_run0', 'r')
    lines = file.readlines()
    # Edit for loop for dfiles
    lines[9] = 'for i in {{{}..{}}}\n'.format(
        i*files_per_script, (i+1)*files_per_script-1)
    if i == scripts:
        lines[9] = 'for i in {{{}..{}}}\n'.format(
            i*files_per_script, i*files_per_script+files_last_scripts)
    file.close()
    filename = 'mis_run{}'.format(i)
    with open(filename, 'w') as f:
        f.write("".join(lines))

    # Run script
    subprocess.call('sbatch {}'.format(filename), shell=True)
