#########################
# sensitivity.py
#########################

'''Creating a graph using matplotlib
to look at the objective flux at different
cut offs of total abundance'''

import pandas as pd
import matplotlib.pyplot as plt 


results_df = pd.read_csv("sensitivity.csv", sep = ",")
print(results_df)

fig, ax = plt.subplots(layout = 'constrained')
ax.plot(results_df['Total Abundance'], results_df['Biomass Flux'], color = 'green', marker = "x")
ax.set_xlabel("Total abundance", fontsize = 18)
ax.set_ylabel("Community Biomass Flux", fontsize = 16)
ax.grid()
plt.show()