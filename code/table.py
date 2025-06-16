#%% Table of φ and α values
'''
Figure 4 (e-f): Table of φ and α values
This script creates a table of φ and α values with color coding, the values are obtained from blocking events data
'''
import matplotlib.pyplot as plt
import numpy as np
import os

def plt_colorbar(data, AX, cmap, label):
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Make your axes dividable
    divider = make_axes_locatable(AX)

    # Create a new axes for colorbar
    cax = divider.append_axes("bottom", size="5%", pad=0.2)

    # Assuming phi_data contains your numerical values
    vmin = np.min(data)
    vmax = np.max(data)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(cmap)  # or whatever colormap you used for phi_colors

    # Create colorbar
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])  # This is required for ScalarMappable
    cbar = fig.colorbar(mappable, orientation='horizontal', cax = cax, shrink = 1, aspect=40, alpha = alph)
    # cbar.set_label(f'{label}', rotation=0, labelpad=0, loc='left')
    cbar.set_label('')  # Remove the default label
    cbar.ax.text(-0.05, 0.5, f'{label}', transform=cbar.ax.transAxes, 
                rotation=0, ha='left', va='center', fontsize = 14)
    return

def font(table):
    for (i, j), cell in table.get_celld().items():
        cell.set_text_props(ha='center', va='center', fontsize=14)
    return

#%% main code
##  Data labels
columns = [ "ERA5", "CESM1 – CTRL", "CESM1 – RCP 8.5"]
rows = [
    "Northern Pacific (JJA)", "Northern Pacific (DJF)",
    "Northern Atlantic (JJA)", "Northern Atlantic (DJF)",
    "Southern Pacific (JJA)", "Southern Pacific (DJF)"
]

# φ values

phi_data = np.array([
    [ 0.71, 0.79, 0.82],
    [ 0.63, 0.78, 0.76],
    [0.66, 0.72, 0.72],
    [ 0.78, 0.83, 0.81],
    [ 0.66, 0.73, 0.68],
    [ 0.61, 0.67, 0.65]
])

# α values
alpha_data = np.array([
    [0.051, 0.041, 0.031],
    [0.054, 0.042, 0.037],
    [0.044, 0.033, 0.035],
    [0.041, 0.030, 0.028],
    [0.031, 0.035, 0.035],
    [0.043, 0.033, 0.026]
])

# Normalize for coloring
phi_norm = (phi_data - phi_data.min()) / (phi_data.max() - phi_data.min())
alpha_norm = (alpha_data - alpha_data.min()) / (alpha_data.max() - alpha_data.min())

# Color maps
phi_colors = plt.cm.Blues(phi_norm)
alpha_colors = plt.cm.Reds(alpha_norm)  # Invert so lower alpha = darker

# Set alpha value for transparency (e.g., 0.5 = 50% transparent)
alph = 0.6
phi_colors[..., -1] = alph
alpha_colors[..., -1] = alph

#%% Create figure
basepath = os.path.expanduser("~/Github")
savingpath = f'{basepath}/plots/Fig4'
os.makedirs(savingpath, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 4))
ax1.axis("off")
ax2.axis("off")

# φ table
phi_table = ax1.table(cellText=np.round(phi_data, 2),
                      rowLabels=rows,
                      colLabels=columns,
                      cellColours=phi_colors,
                      loc='center')

font(phi_table)
phi_table.scale(1, 2.5) 
ax1.set_title("Temporal Correlation, φ", fontsize=14)
plt_colorbar(phi_data, ax1, 'Blues', 'φ')

# α table
alpha_table = ax2.table(cellText=np.round(alpha_data, 3),
                        rowLabels=rows,
                        colLabels=columns,
                        cellColours=alpha_colors,
                        loc='center')
font(alpha_table)
alpha_table.scale(1, 2.5) 
ax2.set_title("Onset Probability, α", fontsize=14)
plt_colorbar(alpha_data, ax2, 'Reds', 'α')

plt.tight_layout()
plt.savefig(f"{savingpath}/Table.png", dpi = 600)
