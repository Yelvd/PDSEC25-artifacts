import sys
import io
import os
import pandas as pd
import numpy as np
import yaml
import sys
import re
from . import misc
from os import listdir
from os.path import isfile, join
from os import walk


def load_archer(experiment_dir, meta):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    results = {}
    for f in files:
        if re.match("energy.csv", f):
            energy_filename = f

    for f in files:
        if ".out" in f and "bull" not in f and "post" not in f:
            out_filename = f

    energy_string = ""
    with open(f"{experiment_dir}/{energy_filename}") as f:
        for l in f:
            if len(l.strip()):
                last_line = l
        last_line = re.sub(' +', ' ', last_line.strip()).split(" ")

        t = last_line[21].split(':')
        results["time"] = float(float(t[0]) * 3600 + float(t[1]) * 60 + float(t[2]))
        results["energy"] = float(last_line[28][:-1])*1000
        results["power"] = results["energy"] / results["time"]

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])

    with open(f"{experiment_dir}/{out_filename}") as f:
        lines = [line.rstrip() for line in f]

        for l in lines:
            if "HemoCell:" in l:
                results["hemocell"] = float(l.split(" ")[-1])
                results["hemocell-energy"] = results["hemocell"] * results["power"]
            if "iterate:" in l:
                results["iterate"] = float(l.split(" ")[-1])
                results["iterate-energy"] = results["iterate"] * results["power"]

    for k in results.keys():
        results[k] = [results[k]]

    df = pd.DataFrame.from_dict(results)
    return df


def load_snellius(experiment_dir, meta):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    results = {}
    energy_filename = None
    energy_2_filename = None
    out_filename = None
    for f in files:
        if re.match("ear-db.csv.*time.csv", f):
            energy_filename = f

        if re.match("ear.csv", f):
            energy_2_filename = f

    for f in files:
        if ".out" in f and "bull" not in f and "post" not in f:
            out_filename = f

    if not (energy_filename or energy_2_filename) or not out_filename:
        print(f"{experiment_dir}: MISSING ear.csv or out file, skipping ....")
        return

    if energy_filename:
        df = pd.read_csv(f"{experiment_dir}/{energy_filename}", delimiter=';')
        results["power"] = df['DC_NODE_POWER_W'].values[0]
        results["time"] = df['TIME_SEC'].values[0]
        results["energy"] = results['power'] * results['time']
    else:
        energy_string = ""
        with open(f"{experiment_dir}/{energy_2_filename}") as f:
            for l in f:
                if len(l.strip()):
                    last_line = l
            last_line = re.sub(' +', ' ', last_line.strip()).split(" ")
            results["power"] = float(last_line[7])
            results["time"] = float(last_line[6])
            results["energy"] = float(last_line[10])

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        print(experiment_dir)
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])

    with open(f"{experiment_dir}/{out_filename}") as f:
        lines = [line.rstrip() for line in f]

        for l in lines:
            if "HemoCell:" in l:
                results["hemocell"] = float(l.split(" ")[-1])
                results["hemocell-energy"] = results["hemocell"] * results["power"]
            if "iterate:" in l:
                results["iterate"] = float(l.split(" ")[-1])
                results["iterate-energy"] = results["iterate"] * results["power"]

    for k in results.keys():
        results[k] = [results[k]]

    df = pd.DataFrame.from_dict(results)
    return df


def load_eviden(experiment_dir, meta):
    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    results = {}
    for f in files:
        if re.match("report.*.csv", f):
            energy_filename = f

    for f in files:
        if ".out" in f and "bull" not in f:
            out_filename = f

    df = pd.read_csv(f"{experiment_dir}/{energy_filename}")
    results["power"] = np.mean(df['Power (W)'])
    results["time"] = np.max(df['Timestamp (s)'])
    results["energy"] = results["power"] * results["time"]

    if meta["General"].get("Benchmark") is None:
        print("WARNING FIXING TO FIXED")
        results["benchmark"] = "fixed"
    else:
        results["benchmark"] = meta["General"]["Benchmark"]

    results["platform"] = meta["General"]["Platform"]

    results["tasks"] = int(meta["General"]["Tasks"])

    with open(f"{experiment_dir}/{out_filename}") as f:
        lines = [line.rstrip() for line in f]

        for l in lines:
            if "HemoCell:" in l:
                results["hemocell"] = float(l.split(" ")[-1])
                results["hemocell-energy"] = results["hemocell"] * results["power"]
            if "iterate:" in l:
                results["iterate"] = float(l.split(" ")[-1])
                results["iterate-energy"] = results["iterate"] * results["power"]

    for k in results.keys():
        results[k] = [results[k]]

    df = pd.DataFrame.from_dict(results)
    return df


def load_das6(experiment_dir, meta):

    print(f"loading {experiment_dir}")

    files = [f for f in listdir(experiment_dir) if isfile(join(experiment_dir, f))]
    results = {}

    for f in files:
        if ".out" in f and "bull" not in f and "post" not in f:
            out_filename = f

    results["benchmark"] = meta["General"]["Benchmark"]
    results["platform"] = meta["General"]["Platform"]
    results["tasks"] = int(meta["General"]["Tasks"])
    results["jobid"] = int(meta["General"]["Id"])

    with open(f"{experiment_dir}/{out_filename}") as f:
        lines = [line.rstrip() for line in f]

        for l in lines:
            if "HemoCell:" in l:
                results["hemocell"] = float(l.split(" ")[-1])
                # results["hemocell-energy"] = results["hemocell"] * results["power"]
            if "iterate:" in l:
                results["iterate"] = float(l.split(" ")[-1])
                # results["iterate-energy"] = results["iterate"] * results["power"]

    nodes = {}
    nodes_results = {}
    for f in files:
        if not "likwid" in f or not "txt" in f:
            continue

        if "tmp" in f:
            continue

        node = f.split('.')[0].split('_')[2]
        task = f.split('.')[0].split('_')[3]
        
        if not node in nodes:
            nodes[node] = []

            with open(f"{experiment_dir}/{f}") as f:
                lines = [line.rstrip() for line in f]
                nodes_results[node] = (float(lines[-1].split(',')[1]),  float(lines[-2].split(',')[1]))

        nodes[node].append(task)

    nnodes = len(nodes.keys())
    for k in results.keys():
        results[k] = [results[k] for _ in range(nnodes)]

    tmp = []
    tmp_e = []
    tmp_p = []
    for k in nodes.keys():
        tmp.append(nodes[k])
        tmp_p.append(nodes_results[k][0])
        tmp_e.append(nodes_results[k][1])

    results["nodes"] = nodes.keys() 
    results["tasks"] = tmp
    results["energy"] = tmp_e
    results["power"] = tmp_p

    # print(results)
    return pd.DataFrame.from_dict(results)

    
def load_experiment(experiment_dir):
    
    with open(f"{experiment_dir}/meta.yaml") as f:
        meta = yaml.load(f, Loader=yaml.FullLoader)

    # if meta["General"]["Platform"] in ["snellius_rome", "snellius_genoa", "snellius-genoa", "snellius-rome"]:
        # return load_snellius(experiment_dir, meta)

    # elif meta["General"]["Platform"] in ["archer_rome", "archer"]:
        # return load_archer(experiment_dir, meta)

    if meta["General"]["Platform"] in ["das6"]:
        return load_das6(experiment_dir, meta)

    else:
        print("Platform not found")
        return None
        


def load_energy(results_dir, experiment=None):
    experiment = misc.translate_exp_to_internal[experiment]
    data = []
    results_dir = f"{results_dir}/"
    w = walk(results_dir)
    for (dirpath, dirnames, filenames) in w:
        if "meta.yaml" not in filenames:
            continue

        with open(f"{dirpath}/meta.yaml") as f:
            meta = yaml.load(f, Loader=yaml.FullLoader)

        data.append(load_experiment(dirpath))

    df = pd.concat(data)

    match experiment:
        case "fixed":
            df = df.loc[df["benchmark"].isin(["fixed", "cube-fixed"])]
        case "fixed-large":
            df = df.loc[df["benchmark"].isin(["fixed-large", "cube-fixed-large"])]
        case "relative":
            df = df.loc[df["benchmark"].isin(["relative", "cube-relative"])]

    # df["benchmark"] = [misc.translate_exp[b] for b in df["benchmark"]]
    # df['platform'] = [misc.translate_platform[p] for p in df['platform']]
    return df 
