from os import listdir
from os.path import isfile, join
from scripts.tools import misc

import sys

figsDir = "artifacts"
scriptDir = "scripts"
matplotlibStyle_ = "matplotlib-style.rc"
resultsDir_ = "results/"
CASES = ["imbalance-particle", "imbalance-fluid"]
SCOREP_EXPS = [f"{exp.path}" for exp in os.scandir(resultsDir_)]

ALL_CASES = ["C1", "C2", "C3", "C3-2"]

PDSEC25_artifacts = [f"{figsDir}/{C}.png" for C in ALL_CASES]


rule PDSEC25:   
    input:
        PDSEC25_artifacts

rule parse_scorep:
    input:
        [f"{exp}/results.csv" for exp in SCOREP_EXPS],


rule parse_scorep_exp:
    input:
        "{expdir}/SCOREP/profile.cubex",
    output:
        "{expdir}/results.csv",
    shell:
        "python3 scripts/exp-to-csv.py -s {wildcards.expdir}/"
        

rule:
    input:
        "scripts/plot-case.py"
    params:
        resultsDir = f"{resultsDir_}/" + "{case}",
        matplotlibStyle = matplotlibStyle_,
        platform = "das6",
        caseName = "{case}"
    output:
        pdf = "{figsDir}/{case}.pdf",
        png = "{figsDir}/{case}.png"
    script:
        "scripts/plot-case.py"
