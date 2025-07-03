import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--cif", dest="cif",
                  help="cif file name, optional")
parser.add_option("--Ecut", dest="Ecut", type=float, default=100.0,
                  help="Energy cutoff, optional")
parser.add_option("--Ncut", dest="Ncut", type=int, default=-1,
                  help="Number cutoff, optional")
(options, args) = parser.parse_args()

# Read and parse the log file
os.system(f"grep data {options.cif} > log.txt")
data = []
with open('log.txt', 'r') as f:
    for line in f:
        if line.startswith('data_'):
            # Split by '-' and extract values
            parts = line.strip().split('-')
            density = float(parts[3][1:])  # Remove 'd' and convert to float
            spg = int(parts[4][3:])        # Remove 'spg' and convert to int
            energy = float(parts[5][1:])*-96.485    # Remove 'e' and convert to float
            data.append({'Density_g_m3': density, 'Space_group': spg, 'Energy_kJ_mol': energy})

# Create DataFrame from parsed data
print("Total number of structures:", len(data))
df = pd.DataFrame(data)

#df = pd.read_csv('crystal_data.csv')
df['Energy_kJ_mol'] -= df['Energy_kJ_mol'].min()  # Normalize energy for better visualization
ids = np.argsort(df['Energy_kJ_mol'])
if options.Ncut > -1: ids = ids[:options.Ncut]

for key in ['Density_g_m3', 'Space_group', 'Energy_kJ_mol']:
    df[key] = df[key][ids]

# Create a figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8))

# Scatter plot in first subplot
scatter = ax1.scatter(df['Density_g_m3'], df['Energy_kJ_mol'], c=df['Space_group'], alpha=0.5, cmap='plasma')
ax1.set_xlabel('Density (g/m$^3$)')
ax1.set_ylabel('Energy (kJ/mol)')
ax1.set_ylim(-3, options.Ecut)
ax1.set_xlim(df['Density_g_m3'].min()-0.05, df['Density_g_m3'].max()+0.05)  # Adjust x-axis limit for better visibility
ax1.set_title('Energy vs Density (colored by Space Group)')
ax1.grid(True)
plt.colorbar(scatter, ax=ax1, label='Space Group', fraction=0.046, pad=0.04)

# Histogram in second subplot
# Group by space group and calculate mean energy and count
sg_stats = df.groupby('Space_group').agg({'Energy_kJ_mol': 'mean', 'Space_group': 'count'}).rename(columns={'Space_group': 'count'})

# Get the actual space group numbers that exist in the data
existing_sg = sorted(df['Space_group'].unique())

# Create mapping of actual space groups to sequential indices
sg_mapping = {sg: i+1 for i, sg in enumerate(existing_sg)}

# Create bar plot with sequential indices
bars = ax2.bar(range(1, len(existing_sg) + 1), sg_stats.loc[existing_sg, 'count'], edgecolor='black', alpha=0.5)

# Set the x-ticks to show the actual space group numbers
ax2.set_xticks(range(1, len(existing_sg) + 1))
ax2.set_xticklabels(existing_sg)

# Create colormap based on mean energy values
norm = plt.Normalize(sg_stats['Energy_kJ_mol'].min(), sg_stats['Energy_kJ_mol'].max())
colors = plt.cm.viridis(norm(sg_stats.loc[existing_sg, 'Energy_kJ_mol']))

# Create bar plot with colors based on mean energy
bars = ax2.bar(range(1, len(existing_sg) + 1), sg_stats.loc[existing_sg, 'count'],
               color=colors, edgecolor='black', alpha=0.7)

# Add colorbar for energy values
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
plt.colorbar(sm, ax=ax2, label='Mean Energy (kJ/mol)', fraction=0.046, pad=0.04)

# Create bar plot with colors based on mean energy
bars = ax2.bar(range(1, len(existing_sg) + 1), sg_stats.loc[existing_sg, 'count'],
               color=colors, edgecolor='black', alpha=0.7)

# Add count numbers and lowest energy on top of each bar
for bar in bars:
    height = bar.get_height()
    sg_num = existing_sg[int(bar.get_x())]  # Get space group number
    min_energy = df[df['Space_group'] == sg_num]['Energy_kJ_mol'].min()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
             #f'{int(height)}\n{min_energy:.1f}',
             f'{min_energy:.1f}',
             ha='center', va='bottom')

ax2.set_xlabel('Space Group Number')
ax2.set_ylabel('Frequency (log scale)')
ax2.set_title(f'Space Group Distribution (Top {len(df)} structures)')
ax2.set_yscale('log')  # Set y-axis to logarithmic scale
ax2.set_ylim(1, df['Space_group'].value_counts().max() * 2.5)  # Adjust y-axis limit for better visibility
ax2.set_xlim(0.5, len(existing_sg) + 0.5)  # Adjust x-axis limits to center bars
plt.tight_layout()
plt.savefig('crystal_data_plot.png', dpi=300)
plt.close()
