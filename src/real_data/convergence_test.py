import numpy as np
import pandas as pd
import pickle
import subprocess
import matplotlib.pyplot as plt
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--static', action='store_true',
                    help='If a static image should be shown instead of the live animation.')
args_params = parser.parse_args()

root = subprocess.check_output(
    'pwd', shell=True).decode('utf-8').replace('\n', '')

# Get number of planets
nplanets = int(root[-12])

# Open parameter names list
parnames = pickle.load(open('parnames.p', 'rb'))
parnames.append('logL')  # Append logL value

# Import data
samples = np.loadtxt(root+f'/hd40307_k{nplanets}_phys_live.txt')
data = pd.DataFrame(samples, columns=parnames)  # Pandas Data Frame

if args_params.static:
    # Static
    plt.figure(2)
    plt.plot(data['planet1_period'], '.')
    plt.title('Planet 1 Period')
    plt.show()
else:
    # Live animation
    plt.ion()
    fig, ax = plt.subplots(nplanets, 1, figsize=(6, 18))
    ax = np.atleast_1d(ax)
    for i in range(120):
        samples = np.loadtxt(root+f'/hd40307_k{nplanets}_phys_live.txt')
        data = pd.DataFrame(samples, columns=parnames)  # Pandas Data Frame
        for n in range(1, nplanets+1):
            ax[n-1].cla()
            ax[n-1].plot(data[f'planet{n}_period'], '.')
            ax[n-1].set_title(
                f'Planet {n} Period, Median = {np.median(data[f"planet{n}_period"]):4.2f}')
            #plt.ylim([-200, 10200])
        plt.tight_layout()
        plt.show()
        plt.pause(10)
