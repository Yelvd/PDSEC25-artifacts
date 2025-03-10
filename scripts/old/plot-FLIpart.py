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
results_dir = f"{snakemake.params['resultsDir']}/das6-imbalance-part"

df_energy = energy.load_energy(results_dir)
# df_energy['freq'] = [float(x.split('-')[1]) for x in df_energy['benchmark']]

df_energy = df_energy.groupby(["jobid"]).agg({"energy": ['sum']})
df_energy = df_energy.reset_index()
df_energy.columns = df_energy.columns.get_level_values(0)


df = scorep.load_scorep(results_dir, region_filter=['iterate'])
df = df.loc[df['platform'] == snakemake.params['platform']]


# change to == 0 with new data
df = df.loc[df['FLIfluid'] == 0]
df = df.loc[df['particles'] > 0]

df = pd.merge(df, df_energy, on="jobid", how="left")

figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[0])
fig, (ax_t, ax_fli, ax_p) = plt.subplots(3, figsize=figsize, sharex=True)

xs = np.unique(df['FLIpart'])
FLIpartes = [[],[]]
Energy = [[],[]]
exs = [[],[]]
comps = [[],[]]
comps_avg = [[],[]]
first = True
for ff in xs:
    tmp_df = df.loc[df['FLIpart'] == ff]
    tmp_df = tmp_df.groupby(["jobid"]).agg({"FLIpart" : ['mean'], "execution": ['mean', 'std', 'max'], "comp": ['mean', 'std', 'max'], "energy": ['mean', 'std']})
    tmp_fli = tmp_df["comp"]['max'] / tmp_df["comp"]['mean'] - 1
    print(tmp_df)

    if first:
        ax_t.scatter([ff] * tmp_df['comp']['max'].size, tmp_df['comp']['max'], color=color_cycle[1], label="comp max")
        ax_t.scatter([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], color=color_cycle[2], label="comp mean")
        first = False
    else:
        ax_t.scatter([ff] * tmp_df['comp']['max'].size, tmp_df['comp']['max'], color=color_cycle[1])
        ax_t.scatter([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], color=color_cycle[2])
    ax_t.errorbar([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], yerr=tmp_df['comp']['std'], color=color_cycle[2])

    FLIpartes[0].append(np.mean(tmp_fli))
    FLIpartes[1].append(np.std(tmp_fli))
    Energy[0].append(np.mean(tmp_df["energy"]["mean"]) / 1000)
    Energy[1].append(np.std(tmp_df["energy"]["mean"]) / 1000)
    exs[0].append(np.mean(tmp_df['execution']['max']))
    exs[1].append(np.std(tmp_df['execution']['max']))
    comps[0].append(np.mean(tmp_df['comp']['max']))
    comps[1].append(np.std(tmp_df['comp']['max']))
    comps_avg[0].append(np.mean(tmp_df['comp']['mean']))
    comps_avg[1].append(np.std(tmp_df['comp']['mean']))
    
    
ax_t.errorbar(xs, exs[0], exs[1], color=color_cycle[0], label="total")
ax_fli.errorbar(xs, FLIpartes[0], FLIpartes[1], color=color_cycle[3], label="FLI")
ax_p.errorbar(xs, Energy[0], Energy[1], color=color_cycle[3], label="FLI")
        # ax_p.scatter([ff] * tmp_df['energy']['mean'].size, tmp_df['energy']['mean'], color=color_cycle[4], label="comp mean")

# ax_t.errorbar(xs, FLIpartes[0], FLIpartes[1], color=color_cycle[1], label="comp max")
# ax_t.errorbar(xs, FLIpartes[0], FLIpartes[1], color=color_cycle[2], label="comp mean")

ax_t.legend()
ax_fli.legend()

ax_t.set_title("9 Node cube with changing particle imbalance DAS6")
ax_p.set_xlabel("Theoretical Fractional Particle Imbalance")
ax_t.set_ylabel("Time [s]")
ax_fli.set_ylabel("FLI")
ax_p.set_ylabel("Energy [kJ]")

ax_t.set_ylim(0)
ax_p.set_ylim(0)
ax_fli.set_ylim(0, 1.5)

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
