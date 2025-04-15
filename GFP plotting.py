import numpy as np #used for numerical operations 
import matplotlib.pyplot as plt #for plotting/visualization of graph
import pandas as pd # for working with strucutred data (data frames)


# Load your dataset (adjust the file path accordingly)
file_path = "/Users/madelinelucas/Desktop/Abreu lab/GFP fitness"
df = pd.read_csv(f'{file_path}/GFP Relative abundance.csv')
print(df.columns)

#explore unique dilutions values in the dataset
print(df['Dilution'].unique())
# Assuming your data is stored in a pandas DataFrame named 'df'
# Example: df = pd.read_csv('your_data.csv')
#%%
#creating a subset of the data containing only the starting point(cycle 0) where dilution =0
df_start=df[df["Dilution"]==0]
# Initial lists to store parsed data from cycle 0
start_sample=[]
start_reps=[]
start_GFP=[]
start_nonGFP=[]

#create a list representing the cylce (alwasy 0 here) repeated for each sample
cycle=[0]*45 #3 samples x 15 wells=45 rows
# build a list of dilution values for cycle 0 manually (5 replicates per sample * 3 dilutions)
start_dilution=[.1]*5
start_dilution.extend([.01]*5)
start_dilution.extend([.001]*5)
start_dilution*=3  #repeate for 3 samples

#dictionary to map each sample (pair) to corresponding 15 replicate well lables
repdict={}
repdict['P1']=['A1','B1','C1','D1','E1','A2','B2','C2','D2','E2','A3','B3','C3','D3','E3']
repdict['P2']=['A4','B4','C4','D4','E4','A5','B5','C5','D5','E5','A6','B6','C6','D6','E6']
repdict['P3']=['A7','B7','C7','D7','E7','A8','B8','C8','D8','E8','A9','B9','C9','D9','E9']

#loop over each sample (P1,P2, P3) and extract data from df_start
for pair in df_start['sample'].unique().tolist():
    this_pair_df=df_start[df_start['sample']==pair]
    start_sample.extend([pair]*15) #repeate the sample name 15 times (P1 5 replicates*3 dilutions)
    start_nonGFP.extend(df_start[df_start['sample']==pair]['GFP - Count/uL'].tolist()*3) #each count is repeated 3 times- this is for 0 growth cycle since no dilutions 
    start_GFP.extend(df_start[df_start['sample']==pair]['GFP + Count/uL'].tolist()*3)
    start_reps.extend(repdict[pair]) #add replicate names
    
#Creates a new DataFrame for cycle 0 using the list above
df_startfull=pd.DataFrame({
    'sample': start_sample,
    'Cycle':cycle,
    'Dilution':start_dilution,
    'Replicates': start_reps,
    'GFP - Count/uL': start_nonGFP,
    'GFP + Count/uL':start_GFP
    
})

#Calculates GFP fraction for the start (cycle 0) dataset
df_startfull['GFPfrac']=df_startfull['GFP + Count/uL']/(df_startfull['GFP + Count/uL']+df_startfull['GFP - Count/uL'])

#Also calculates GFp fraction in the full dataset
df['GFPfrac']=df['GFP + Count/uL']/(df['GFP + Count/uL']+df['GFP - Count/uL'])
#Set 'sample' column as the index for easier grouping later
df=df.set_index('sample')
df_startfull=df_startfull.set_index('sample')

#%%
# Creating a subset of df excluding the zero dilution entries (which we already handled)  
#COLON AFTER THE FOR FOR looping
sub_df=df[['Cycle','Dilution','Replicates','GFP - Count/uL','GFP + Count/uL','GFPfrac']]
sub_df=sub_df[sub_df['Dilution']!=0] #only take data from cycle >0

#combine this subset with the prepared cycle 0 data
combined_df=pd.concat([sub_df, df_startfull])
#%%
#creates a 3x3 grid of subplots: 3 samples x 3 dilution 
fig, axs = plt.subplots(3, 3, figsize=(12, 9), sharex=True, sharey=True)

#define the sample groups and dilution levels
pairs = ['P1', 'P2', 'P3']
dils = df['Dilution'].unique().tolist()[1:]#skip dilution 0 (already handled)

#loop through each sample (row of subplot grid)
for k in range(3):
    pair = pairs[k]
    pair_df = combined_df[combined_df.index == pair] #filter data for this sample
    
#loop through each dilution (columbn of subplot grid)
    for j in range(3):
        dil = dils[j]
        dil_df = pair_df[pair_df['Dilution'] == dil] #filter for this dilution
        
        ax = axs[k, j] #select appropriate subplot
        
    #group by replicate and plot GFP fraction over time (cycle)
        
        for label, group in dil_df.groupby('Replicates'):
            sorted_group = group.sort_values('Cycle')
            ax.plot(sorted_group['Cycle'], sorted_group['GFPfrac'], marker='o', label=f"Rep {label}")
        
        ax.set_ylim(0, 1) # Kepp y-axis consistent fro all plots

        # Only set column titles on the top row a
        if k == 0:
            ax.set_title(f'Dilution {dil}')
        
        # Only set y-axis labels on the first column
        if j == 0:
            ax.set_ylabel(f'Pair {pair}')

# Optionally, add a common X and Y label
fig.supxlabel('Cycles')
fig.supylabel('GFP Fraction of pairs')

# Add legend only once (optional)
# You can also use fig.legend() if you want a shared one
# axs[0, 0].legend(loc='upper right')

plt.tight_layout()
plt.show()

            
            