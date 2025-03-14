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
results_dir = snakemake.params['resultsDir']

df_energy = energy.load_energy(results_dir)

runs=['0', '1', '2', '3', '4', '5']


# df_energy = df_energy.groupby(["jobid"]).agg({"energy": ['sum']})
# df_energy = df_energy.reset_index()
# df_energy.columns = df_energy.columns.get_level_values(0)


df = scorep.load_scorep(results_dir)

fluid_regions = ['collideandstream', 'setexternalvector']
part_regions = ['advanceparticles', 'applyconstitutivemodel', 'deletenonlocalparticles', 'syncenvelopes']
coupl_regions = ['interpolatefluidvelocity', 'spreadparticleforce']

# new_rows = []
# for i in np.unique(df['jobid']):
    # tmp_df = df.loc[df['jobid'] == i]
    # for t in np.unique(tmp_df[' Thread ID']):
        # tmp_df_2 = tmp_df.loc[tmp_df[' Thread ID'] == t]

        # r = tmp_df_2.reset_index().loc[0].copy()
        # tmp = tmp_df_2.loc[tmp_df_2['regionName'].isin(fluid_regions)]
        # r['regionName'] = 'fluid_component'
        # r['mpi'] = np.sum(tmp['mpi'])
        # r['comp'] = np.sum(tmp['comp'])
        # r['execution'] = np.sum(tmp['execution'])
        # new_rows.append(r)

        # r = tmp_df_2.reset_index().loc[0].copy()
        # tmp = tmp_df_2.loc[tmp_df_2['regionName'].isin(part_regions)]
        # r['regionName'] = 'part_component'
        # r['mpi'] = np.sum(tmp['mpi'])
        # r['comp'] = np.sum(tmp['comp'])
        # r['execution'] = np.sum(tmp['execution'])
        # new_rows.append(r)

        # r = tmp_df_2.reset_index().loc[0].copy()
        # tmp = tmp_df_2.loc[tmp_df_2['regionName'].isin(coupl_regions)]
        # r['regionName'] = 'couple_component'
        # r['mpi'] = np.sum(tmp['mpi'])
        # r['comp'] = np.sum(tmp['comp'])
        # r['execution'] = np.sum(tmp['execution'])
        # new_rows.append(r)

        
df = df.loc[df['regionName'] == 'iterate']
# df_comps = pd.DataFrame(new_rows)


# df = pd.merge(df, df_energy, on="jobid", how="left")

figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[0])
fig, (ax_t, ax_fli, ax_p) = plt.subplots(3, figsize=figsize, sharex=True)

xs = [0,1,2,3,4,5]
FLIpartes = [[],[]]
Energy = [[],[]]
exs = [[],[]]
comps = [[],[]]
first = True

for ff in xs:
    ids = np.unique(df_energy.loc[df_energy['benchmark'].str.contains(runs[ff])]['jobid'])
    tmp_df = df.loc[df['jobid'].isin(ids)]
    tmp_df = tmp_df.groupby(["jobid"]).agg({"FLIpart" : ['mean'], "execution": ['mean', 'std', 'max'], "comp": ['mean', 'std', 'max']})
    tmp_df = tmp_df.reset_index()
    tmp_fli = tmp_df["comp"]['max'] / tmp_df["comp"]['mean'] - 1

    if first:
        ax_t.scatter([ff] * tmp_df['comp']['max'].size, tmp_df['comp']['max'], color=color_cycle[1], label="comp max")
        ax_t.scatter([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], color=color_cycle[2], label="comp mean")
        first = False
    else:
        ax_t.scatter([ff] * tmp_df['comp']['max'].size, tmp_df['comp']['max'], color=color_cycle[1])
        ax_t.scatter([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], color=color_cycle[2])
    ax_t.errorbar([ff] * tmp_df['comp']['mean'].size, tmp_df['comp']['mean'], yerr=tmp_df['comp']['std'], color=color_cycle[2])


    power = []
    energy = []
    for i in ids:
        power.append(np.sum(df_energy.loc[df_energy['jobid'] == i]['power']))
        energy.append(power[-1] * tmp_df.loc[tmp_df['jobid'] == i]['execution']['max'].to_list()[0])


    Energy[0].append(np.mean(energy) / 1000)
    Energy[1].append((np.min(energy) / 1000) - Energy[0])
    Energy[2].append((np.max(energy) / 1000) - Energy[0])
    FLIpartes[0].append(np.mean(tmp_fli))
    FLIpartes[1].append(np.min(tmp_fli) - FLIpartes[0])
    FLIpartes[2].append(np.max(tmp_fli) - FLIpartes[0])
    exs[0].append(np.mean(tmp_df['execution']['max']))
    exs[1].append(np.min(tmp_df['execution']['max']) - exs[0])
    exs[2].append(np.max(tmp_df['execution']['max']) - exs[0])
    comps[0].append(np.mean(tmp_df['comp']['max']))
    comps[1].append(np.min(tmp_df['comp']['max'])- comps[0])
    comps[2].append(np.max(tmp_df['comp']['max'])- comps[0])
    
    
ax_t.errorbar(xs, exs[0], yerr=[exs[1],exs[2]], color=color_cycle[0], label="total")
ax_fli.errorbar(xs, FLIpartes[0], FLIpartes[1], color=color_cycle[3], label="FLI")
ax_p.errorbar(xs, Energy[0], Energy[1], color=color_cycle[3], label="FLI")
for x,y in zip(xs, Energy[0]):
    plt.annotate(f'{y:.0f}', # this is the text
                     (x,y), # these are the coordinates to position the label
                     textcoords="offset points", # how to position the text
                     xytext=(0,10), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center
            # ax_p.scatter([ff] * tmp_df['energy']['mean'].size, tmp_df['energy']['mean'], color=color_cycle[4], label="comp mean")

# ax_t.errorbar(xs, FLIpartes[0], FLIpartes[1], color=color_cycle[1], label="comp max")
# ax_t.errorbar(xs, FLIpartes[0], FLIpartes[1], color=color_cycle[2], label="comp mean")

ax_t.legend()
ax_fli.legend()

ax_t.set_title("16 Node cube with changing particle imbalance DAS6")
ax_p.set_xlabel("Run")
ax_t.set_ylabel("Time [s]")
ax_fli.set_ylabel("FLI")
ax_p.set_ylabel("Energy [kJ]")

ax_t.set_ylim(0, 300)
ax_p.set_ylim(0, 400)
ax_fli.set_ylim(0, 2.2)
ax_fli.set_xticks(xs, ["B", "IB", "IB-2.8", "S1", "S2", "S3"])

plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
