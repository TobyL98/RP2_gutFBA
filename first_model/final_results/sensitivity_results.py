#########################
# sensitivity_results.py
#########################

'''Creating a graph using matplotlib
to look at the objective flux at different
cut offs of total abundance'''

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path 

def sensitivity_plotter(ax, results_df, color):

    """Function that plots the correct plot"""

    ax.plot(results_df['Total Abundance'], results_df['Biomass Flux'], color = color, marker = "x")
    ax.grid()

# read in data
results_healthy_path = Path("sensitivity_healthy/sens_healthy_results.csv")
results_healthy_df = pd.read_csv(results_healthy_path, sep = ",")

results_CRC_path = Path("sensitivity_CRC/sens_CRC_results.csv")
results_CRC_df = pd.read_csv(results_CRC_path, sep = ",")

# set up subplots and run plots function
fig, (ax1, ax2) = plt.subplots(1, 2, layout = "constrained", sharey = True)
sensitivity_plotter(ax1, results_healthy_df, color = "green")
sensitivity_plotter(ax2, results_CRC_df, color = "red")

ax1.set_ylabel("Community Biomass Flux", fontsize = 16)
ax1.set_title("Healthy", fontsize = 20)

ax2.set_title("CRC", fontsize = 20)

fig.supxlabel("Total Abundance", fontsize = 18)
#fig.text(0.5, 0.02, 'Total Abundance', ha='center', fontsize = 18)

plt.show()
#plt.savefig("sensitivity_fig.png")