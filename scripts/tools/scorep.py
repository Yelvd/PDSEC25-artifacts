import sys
import os
import pandas as pd
import numpy as np
import yaml
import sys
import glob
from bs4 import BeautifulSoup
from . import misc
from os import listdir
from os.path import isfile, join
from os import walk

def load_scorep(results_dir, experiment=None, region_filter=None):
    data = []
    # experiment = misc.translate_exp_to_internal[experiment]
    w = walk(results_dir)
    for (dirpath, dirnames, filenames) in w:
        if "meta.yaml" not in filenames:
            continue

        with open(f"{dirpath}/meta.yaml") as f:
            meta = yaml.load(f, Loader=yaml.FullLoader)

        if meta["General"]["Monitor-tool"] in ["scorep"]:
            data.append(load_experiment(dirpath, region_filter))

    df = pd.concat(data)

    if experiment is not None:
        df = df.loc[df["benchmark"] == experiment]

    return df 


def load_experiment(experiment_dir, region_filter):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    
    results = {}
    for f in files:
        if ".out" in f and "bull" not in f and "post" not in f:
            out_filename = f

    with open(f"{experiment_dir}/{out_filename}") as f:
        lines = [line.rstrip() for line in f]

        for l in lines:
            if "Smallest atomic-block:" in l:
                results["small-atmoic-block"] = l.split(" ")[-1]
            if "Largest atomic-block:" in l:
                results["large-atmoic-block"] = l.split(" ")[-1]
            if "(main)   nCells (global)" in l:
                results["particles"] = int(l.split(" ")[-1])
    
    with open(f"{experiment_dir}/meta.yaml") as f:
        meta = yaml.load(f, Loader=yaml.FullLoader)

    results["benchmark"] = meta["General"]["Benchmark"]
    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["nodes"] = int(meta["General"]["Nodes"])
    results["jobid"] = int(meta["General"]["Id"])


    filename = glob.glob(f'{experiment_dir}/config*.xml')[0]
    # Reading the data inside the xml file to a variable under the name  data
    with open(filename, 'r') as file:
        xml = file.read()

    # Passing the stored data inside the beautifulsoup parser
    bs_data = BeautifulSoup(xml, 'xml')

    # Using find() to extract attributes of the first instance of the tag
    b_name = bs_data.find('FLIfluid')
    for test in b_name.children:
        results['FLIfluid'] = float(test)

    try:
        # Using find() to extract attributes of the first instance of the tag
        b_name = bs_data.find('FLIpart')
        for test in b_name.children:
            results['FLIpart'] = float(test)
    except:
        results['FLIpart'] = -1


    try:
        data = pd.read_csv(f"{experiment_dir}/results.csv")
    except:
        print(f"could not load {experiment_dir}/results.csv")
        return pd.DataFrame([],[])

    if region_filter is not None:
        data = data.loc[data['regionName'].isin(region_filter)]

    try:
        results["benchmark"] = misc.translate_exp[results["benchmark"]]
    except:
        pass

    try:
        results['platform'] = misc.translate_platform[results['platform']]
    except:
        pass

    for k in results.keys():
        data[k] = results[k]
    
    return data
