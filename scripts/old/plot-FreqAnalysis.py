import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import walk

from tools import misc
from tools import scorep
from tools import energy

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = f"{snakemake.params['resultsDir']}/das6-freq-analysis"

df = scorep.load_scorep(results_dir, region_filter=['iterate'])

df = df.loc[df['platform'] == snakemake.params['platform']]

df_energy = energy.load_energy(results_dir)
df_energy['freq'] = [float(x.split('-')[1]) for x in df_energy['benchmark']]
print(df_energy)

df_energy = df_energy.groupby(["jobid","freq"]).agg({"energy": ['sum', 'mean', 'std', 'max', "min"], "power": ['sum', 'mean', 'std', 'max']})
df_energy = df_energy.reset_index()
df_energy = df_energy.sort_values("freq")
print(df_energy)


figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[0])
fig, (ax_t, ax_e, ax_p) = plt.subplots(3, figsize=figsize, sharex=True)

df['freq'] = [float(x.split('-')[1]) for x in df['benchmark']]

xs = np.unique(df['freq'])
FLIfluides = [[],[]]
exs = [[],[]]
comps = [[],[]]
comps_avg = [[],[]]
first = True
barwidth = .1
e_scale = 1000

for ff in xs:
    tmp_df = df.loc[df['freq'] == ff]
    tmp_df = tmp_df.groupby(["jobid"]).agg({"freq" : ['mean'], "execution": ['mean', 'std', 'max'], "comp": ['mean', 'std', 'max']})
    tmp_fli = tmp_df["comp"]['max'] / tmp_df["comp"]['mean'] - 1
    tmp_e = df_energy.loc[df_energy["freq"] == ff]

    if first:
        ax_t.scatter([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], color=color_cycle[2], label="Average Compute cost")
        ax_e.bar([ff], [tmp_e['energy']['sum'].mean() / e_scale], label="Total PKG Energy", color=color_cycle[3], width=barwidth)
        ax_e.errorbar([ff], [tmp_e['energy']['sum'].mean()/ e_scale], tmp_e['energy']['sum'].std()/ e_scale, color='k')
        ax_p.bar([ff], [tmp_e['power']['sum'].mean()], label="Total PKG Power", color=color_cycle[4], width=barwidth)
        ax_p.errorbar([ff], [tmp_e['power']['sum'].mean()], tmp_e['power']['sum'].std(), color='k')
        first = False
    else:
        ax_t.scatter([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], color=color_cycle[2])
        ax_e.bar([ff], [tmp_e['energy']['sum'].mean()/ e_scale], color=color_cycle[3], width=barwidth)
        ax_e.errorbar([ff], [tmp_e['energy']['sum'].mean()/ e_scale], tmp_e['energy']['sum'].std()/ e_scale, color='k')
        ax_p.bar([ff], [tmp_e['power']['sum'].mean()], color=color_cycle[4], width=barwidth)
        ax_p.errorbar([ff], [tmp_e['power']['sum'].mean()], tmp_e['power']['sum'].std(), color='k')
    ax_t.errorbar([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], yerr=tmp_df['comp']['std'], color=color_cycle[2])

    FLIfluides[0].append(np.mean(tmp_fli))
    FLIfluides[1].append(np.std(tmp_fli))
    exs[0].append(np.mean(tmp_df['execution']['max']))
    exs[1].append(np.std(tmp_df['execution']['max']))
    comps[0].append(np.mean(tmp_df['comp']['max']))
    comps[1].append(np.std(tmp_df['comp']['max']))
    comps_avg[0].append(np.mean(tmp_df['comp']['mean']))
    comps_avg[1].append(np.std(tmp_df['comp']['mean']))
    
    
ax_t.errorbar(xs, exs[0], exs[1], color=color_cycle[0], label="Total runtime")

# ax_t.errorbar(xs, FLIfluides[0], FLIfluides[1], color=color_cycle[1], label="comp max")
# ax_t.errorbar(xs, FLIfluides[0], FLIfluides[1], color=color_cycle[2], label="comp mean")

ax_t.legend()
ax_e.legend()
ax_p.legend()
ax_e.set_ylabel("Energy [kJ]")
ax_p.set_ylabel("Power [W]")

ax_t.set_title("8 Node cube with changing freq DAS6")
ax_t.set_ylabel("Time [s]")

ax_p.set_xlabel("Frequency [GHz]")
ax_p.set_xticks(xs)
ax_t.set_ylim(0)
ax_e.set_ylim(0)
ax_p.set_ylim(0)

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
