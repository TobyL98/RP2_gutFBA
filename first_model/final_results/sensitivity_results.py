#########################
# sensitivity_results.py
#########################

'''Creating a graph using matplotlib
to look at the objective flux at different
cut offs of total abundance'''

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from pathlib import Path 

def sensitivity_plotter(ax, results_df, color):

    """Function that plots the biomass flux plot"""

    ax.plot(results_df['Percentage Abundance'], results_df['Biomass Flux'], color = color, marker = "x")
    ax.grid()

def time_taken_plotter(ax, results_df, color):

    '''Function that plots the time takne plot'''
    ax.plot(results_df['Percentage Abundance'], results_df['Time Taken'], color = color, marker = "x")
    ax.grid()

def label_function(ax, label):

    ax.text(
        0.0, 1.0, label, transform=(
            ax.transAxes + ScaledTranslation(-20/72, +7/72, fig.dpi_scale_trans)),
        fontsize=18, va='bottom')

# read in data
results_healthy_path = Path("sensitivity_healthy/sens_healthy_results.csv")
results_healthy_df = pd.read_csv(results_healthy_path, sep = ",")
results_healthy_df['Percentage Abundance'] = results_healthy_df['Total Abundance'] * 100
print(results_healthy_df.head())

results_CRC_path = Path("sensitivity_CRC/sens_CRC_results.csv")
results_CRC_df = pd.read_csv(results_CRC_path, sep = ",")
results_CRC_df['Percentage Abundance'] = results_healthy_df['Total Abundance'] * 100

# set up subplots and run plots function
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, layout = "constrained")
sensitivity_plotter(ax1, results_healthy_df, color = "green")
sensitivity_plotter(ax2, results_CRC_df, color = "red")
time_taken_plotter(ax3, results_healthy_df, color = "green")
time_taken_plotter(ax4, results_CRC_df, color = "red")

ax1.set_ylabel("Community Biomass Flux\n (mmol/d)", fontsize = 16)
ax3.set_ylabel("Time Taken (mins)", fontsize = 16)
ax1.set_title("Healthy", fontsize = 20)

ax2.set_title("CRC", fontsize = 20)

# set tag labels (A, B, C and D)
label_function(ax1, label= "A")
label_function(ax2, label= 'B')
label_function(ax3, label= "C")
label_function(ax4, label= "D")

fig.supxlabel("Total Abundance (%)", fontsize = 22)
#fig.text(0.5, 0.02, 'Total Abundance', ha='center', fontsize = 18)

plt.show()
#plt.savefig("sensitivity_fig2.png")