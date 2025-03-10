import json
import numpy as np
from os import listdir
from os.path import isfile, join
from os import walk



# Import module
import os

# Assign directory
directory = "results"

fli = []
# Iterate over files in directory
for jobid in os.listdir(directory):
    # Open file
    if not os.path.isfile(f"results/{jobid}/tmp_1/log_1/logfile"):
        continue

    sizes = []
    for stats in os.listdir(f"results/{jobid}/tmp_1/log_1"):
        if "statistics." not in stats:
            continue


        # Open and read the JSON file
        with open(f"results/{jobid}/tmp_1/log_1/{stats}", 'r') as file:
            data = json.load(file)

        # Print the data
        for threadid in data.keys():
            sizes.append(int(data[threadid]['Metrics']['Atomic Block Size']))

    fli.append(np.round(np.max(sizes) / np.mean(sizes) - 1, 4))
    print(f"{jobid} : max = {np.max(sizes)}, avg = {np.mean(sizes)}, fli = {np.round(np.max(sizes) / np.mean(sizes) - 1, 4)}")

fli.sort()
print(fli)
