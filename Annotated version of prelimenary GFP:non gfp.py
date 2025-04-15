import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Only keep rows where Dilution is 0
df_start = df[df["Dilution"] == 0].copy()

# Build a dilution pattern to assign based on known structure
dilutions = np.tile(np.repeat([0.1, 0.01, 0.001], 5), 3)

# Create a mapping from sample to replicate wells
repdict = {
    'P1': ['A1','B1','C1','D1','E1','A2','B2','C2','D2','E2','A3','B3','C3','D3','E3'],
    'P2': ['A4','B4','C4','D4','E4','A5','B5','C5','D5','E5','A6','B6','C6','D6','E6'],
    'P3': ['A7','B7','C7','D7','E7','A8','B8','C8','D8','E8','A9','B9','C9','D9','E9'],
}

# Build the DataFrame
records = []
df_reset = df.reset_index()  # brings 'sample' back as a column
df_start = df_reset[df_reset["Dilution"] == 0].copy()


for pair in df_start['sample'].unique():
    # Extract data for this pair
    pair_data = df_start[df_start['sample'] == pair][['GFP + Count/uL', 'GFP - Count/uL']].values.flatten()
    gfp_counts = np.tile(pair_data[::2], 3)  # +GFP counts repeated 3 times (for 3 cycles)
    nongfp_counts = np.tile(pair_data[1::2], 3)
    reps = repdict[pair]
    cycles = [0] * 45
    
    for i in range(45):
        records.append({
            'sample': pair,
            'Cycle': cycles[i],
            'Dilution': dilutions[i],
            'Replicates': reps[i],
            'GFP - Count/uL': nongfp_counts[i],
            'GFP + Count/uL': gfp_counts[i],
        })

df_startfull = pd.DataFrame(records)
df_startfull['GFPfrac'] = df_startfull['GFP + Count/uL'] / (
    df_startfull['GFP + Count/uL'] + df_startfull['GFP - Count/uL']
)
#%%
from itertools import product

fig, axs = plt.subplots(3, 3, figsize=(12, 9), sharex=True, sharey=True)
pairs = ['P1', 'P2', 'P3']
dils = sorted(df['Dilution'].unique()[1:])

for (k, pair), (j, dil) in product(enumerate(pairs), enumerate(dils)):
    ax = axs[k, j]
    pair_dil_df = combined_df.loc[pair]
    if isinstance(pair_dil_df, pd.DataFrame):  # Handle if only one row
        pair_dil_df = pair_dil_df[pair_dil_df['Dilution'] == dil]
        for label, group in pair_dil_df.groupby('Replicates'):
            ax.plot(group['Cycle'], group['GFPfrac'], marker='o')
    
    if k == 0:
        ax.set_title(f'Dilution {dil}')
    if j == 0:
        ax.set_ylabel(f'Pair {pair}')

fig.supxlabel('Cycle')
fig.supylabel('GFP Fraction')
plt.tight_layout()
plt.show()

