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

def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))

runs=['0', '1', '2', '3', '4', '5']


# df_energy = df_energy.groupby(["jobid"]).agg({"energy": ['sum']})
# df_energy = df_energy.reset_index()
# df_energy.columns = df_energy.columns.get_level_values(0)


df = scorep.load_scorep(results_dir)
df = df.loc[df['regionName'] == 'iterate']

figsize = plt.rcParams['figure.figsize']
figsize = (figsize[0], figsize[1])
fig, (ax_t, ax_fli, ax_p) = plt.subplots(3, figsize=figsize, sharex=True)

xs = [0,1,2,3,4,5]
FLIpartes = [[],[]]
Energy = [[],[]]
exs = [[],[]]
comps = [[],[]]
first = True
size = 3

for ff in xs:
    ids = np.unique(df_energy.loc[df_energy['benchmark'].str.contains(runs[ff])]['jobid'])
    tmp_df = df.loc[df['jobid'].isin(ids)]

    v = ax_t.violinplot(tmp_df['comp'], positions = [ff], showmeans=True, widths=0.8)

    for pc in v['bodies']:
        pc.set_facecolor(color_cycle[1])

    for partname in ('cbars','cmins','cmaxes','cmeans'):
        vp = v[partname]
        vp.set_edgecolor(color_cycle[1])
        vp.set_linewidth(1)
        # pc.set_edgecolor('black')
    # if first:
    # else:
        # ax_t.violinplot(tmp_df['comp'], positions = [ff], showmeans=True, width=0.8)

    tmp_df = tmp_df.groupby(["jobid"]).agg({"FLIpart" : ['mean'], "execution": ['mean', 'std', 'max'], "comp": ['mean', 'std', 'max']})
    tmp_df = tmp_df.reset_index()
    tmp_fli = tmp_df["comp"]['max'] / tmp_df["comp"]['mean'] - 1


    power = []
    energy = []
    for i in ids:
        power.append(np.sum(df_energy.loc[df_energy['jobid'] == i]['power']))
        energy.append(power[-1] * tmp_df.loc[tmp_df['jobid'] == i]['execution']['max'].to_list()[0])


    Energy[0].append(np.mean(energy) / 1000)
    Energy[1].append(np.std(energy) / 1000)
    FLIpartes[0].append(np.mean(tmp_fli))
    FLIpartes[1].append(np.std(tmp_fli))
    exs[0].append(np.mean(tmp_df['execution']['max']))
    exs[1].append(np.std(tmp_df['execution']['max']))
    comps[0].append(np.mean(tmp_df['comp']['max']))
    comps[1].append(np.std(tmp_df['comp']['max']))
    
    
ax_t.errorbar(xs, exs[0], exs[1], color=color_cycle[0], marker='.', fmt='--', lw=1, label="Total")
ax_t.bar(1, -2, color=color_cycle[1], label="Compute")

ax_fli.errorbar(xs, FLIpartes[0], FLIpartes[1], color='k', linestyle='', label="FLI")
ax_fli.bar(xs, FLIpartes[0], color=color_cycle[2], width=0.4, label="$f_{li}$")
for x,y in zip(xs, FLIpartes[0]):
    ax_fli.annotate(f'{y:.1f}', # this is the text
                     (x,y), # these are the coordinates to position the label
                     textcoords="offset points", # how to position the text
                     xytext=(0,4), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center

# ax_fli.errorbar(xs, FLIpartes[0], FLIpartes[1], color='r', fmt='', linestyle='')
ax_p.errorbar(xs, Energy[0], Energy[1], color='k', linestyle='', label="FLI")
ax_p.bar(xs, Energy[0], color=color_cycle[3], width=0.4, label="")
for x,y in zip(xs, Energy[0]):
    plt.annotate(f'{y:.0f}', # this is the text
                     (x,y), # these are the coordinates to position the label
                     textcoords="offset points", # how to position the text
                     xytext=(0,4), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center
            # ax_p.scatter([ff] * tmp_df['energy']['mean'].size, tmp_df['energy']['mean'], color=color_cycle[4], label="comp mean")

# Shrink current axis by 20%
# box = ax_t.get_position()
# ax_t.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
# ax_t.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)

ax_t.legend(ncol=2)
# ax_p.set_xlabel("Run")
ax_t.set_ylabel("Time [s]")
ax_fli.set_ylabel("Imbalance [$f_{li}$]")
ax_p.set_ylabel("Energy [kJ]")

ax_t.set_ylim(0, 320)
ax_p.set_ylim(0, 520)
ax_fli.set_ylim(0, 2.2)
# ax_fli.set_xticks(xs, ["Balanced", "Imbalanced", "Imbalanced-2.8GHz", "S1", "S2", "S3"])
ax_fli.set_xticks(xs, ["B", "IB", "IB@2.8", "S1", "S2", "S3"])
# xt = plt.xticks(rotation=45, ha='right', rotation_mode='anchor' )
# xt.set_in_layout(False) 

fig.align_ylabels((ax_t, ax_fli, ax_p))
plt.savefig(snakemake.output["png"], dpi=300)
plt.savefig(snakemake.output["pdf"])
