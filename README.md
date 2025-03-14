# PDSEC25-artifacts
[![DOI](https://zenodo.org/badge/946292613.svg)](https://doi.org/10.5281/zenodo.15003097)

A repository with the data and script for creating the artifacts for the PDSEC25 paper "Embracing Load Imbalance for Energy Optimizations: a Case-Study"

The results in this repostiry are obtained by running [HemoCell](hemocell.eu) with benchnmarks from [github.com/UvaCsl/Hemocell-Performance-Benchmarks](github.com/UvaCsl/Hemocell-Performance-Benchmarks) on the Distruted ASCII Supercomputer [DAS6](https://www.cs.vu.nl/das/)

## Repository
- `results/`: Contains all the results and metadata for each experiment.
- `scrips/`: A collection of scripts for loading the data and generating the plots.
- `setup/`: A set of scripts used to create and run the experiments.
    

## How to create artifacts

### requirements
1. poetry (https://python-poetry.org/)

### steps

1. `poetry install`
2. `mkdir artifacts`
3. `poetry run snakemake -c 1 PDSEC25`

The plots will be generated in `artifacts/`
