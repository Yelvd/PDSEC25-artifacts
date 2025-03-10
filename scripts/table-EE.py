import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import walk

from tools import misc
from tools import scorep
from tools import energy as energyTool

plt.style.use(snakemake.params["matplotlibStyle"])
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
results_dir = snakemake.params['resultsDir']

runs=['0', '1', '2', '3', '4', '5']
Names=["B", "IB", "IB-2.8", "S1", "S2", "S3"]

fluid_regions = ['collideandstream', 'setexternalvector']
part_regions = ['advanceparticles', 'applyconstitutivemodel', 'deletenonlocalparticles', 'syncenvelopes']
coupl_regions = ['interpolatefluidvelocity', 'spreadparticleforce']
xs = [0,1,2,3,4,5]
iterations = 500

Columns = ['Case'] + Names

data = []

for C in snakemake.params['cases']:
    df_energy = energyTool.load_energy(results_dir + f"/{C}")
    df = scorep.load_scorep(results_dir + f"/{C}")

    df = df.loc[df['regionName'] == 'iterate']

    Energy = [[],[]]

    for ff in xs:
        ids = np.unique(df_energy.loc[df_energy['benchmark'].str.contains(runs[ff])]['jobid'])
        tmp_df = df.loc[df['jobid'].isin(ids)]
        tmp_df = tmp_df.groupby(["jobid"]).agg({"FLIpart" : ['mean'], "execution": ['mean', 'std', 'max'], "comp": ['mean', 'std', 'max']})
        tmp_df = tmp_df.reset_index()
        tmp_fli = tmp_df["comp"]['max'] / tmp_df["comp"]['mean'] - 1

        power = []
        energy = []
        for i in ids:
            power.append(np.sum(df_energy.loc[df_energy['jobid'] == i]['power']))
            energy.append(power[-1] * tmp_df.loc[tmp_df['jobid'] == i]['execution']['max'].to_list()[0])


        Energy[0].append(np.mean(np.array( energy) / iterations))
        Energy[1].append(np.std(np.array( energy ) / iterations))

    data.append([C] + Energy[0])

new_df = pd.DataFrame(data, columns=Columns)
new_df.reset_index().set_index(['Case'])
print(new_df)

i = 0
for c in np.array(new_df.columns):
    if c != "Case":
        new_df[c] = new_df[c].map(lambda x: '{0:.0f}'.format(x))



# print(new_df.to_latex(index=False, hrules=True, column_format='lrrrrrr'))
latex_str = new_df.style.format().hide(axis="index").to_latex(hrules=True, column_format='@{}lrrrrr@{}')

with open(snakemake.output["tex"], "w") as file1:
    file1.write(latex_str)
