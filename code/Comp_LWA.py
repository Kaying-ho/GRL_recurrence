#%% Hovmoller Composite of LWA ##
'''
Figure 2: Composites of Hovmöller diagram for LWA
This script generates composites of Hovmöller diagram for LWA during blocking events in different seasons and regions.
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp
import os

def plot_levels(season):
    '''plotting levels for each season time'''
    if season == "summer":
        lev3 = np.linspace(5e8, 1.7e9, num=13) # Red shading: Current
        lev0 = np.linspace(0.05, 0.29, num=13)  # Purple shading: Next/Previous
    elif season == "winter":
        lev3 = np.linspace(1e9, 3.4e9, num=13)  # Red shading: Current
        lev0 = np.linspace(0.05, 0.31, num=14) # Purple shading: Next/Previous
    else:
        raise ValueError("Invalid season. Choose 'JJA' or 'DJF'.")
    return lev3, lev0

def perform_blockwise_ttest(data, block_size=4, null_hypothesis_mean=0.5, significance_level=0.05):
    num_blocks_y = data.shape[0] // block_size
    num_blocks_x = data.shape[1] // block_size

    p_values = np.empty((num_blocks_y, num_blocks_x))

    for i in range(num_blocks_y):
        for j in range(num_blocks_x):
            block = data[i*block_size:(i+1)*block_size, j*block_size:(j+1)*block_size]
            t_stat, p_val = ttest_1samp(block.reshape(-1), null_hypothesis_mean, alternative='greater')
            p_values[i, j] = p_val

    significant_points = p_values < significance_level
    
    return p_values, significant_points, num_blocks_y, num_blocks_x

def plot_significant_points(significant_points, block_size, nlon, ndays, label=False):
    num_blocks_y, num_blocks_x = significant_points.shape
    for i in range(num_blocks_y):
        for j in range(num_blocks_x):
            if significant_points[i, j]:
                x_center = j * block_size + block_size // 2 - int(nlon / 2)
                y_center = i * block_size + block_size // 2 - int(ndays / 2)
                plt.scatter(
                    x_center, y_center,
                    color='k',
                    marker='o',
                    label='Significant Points' if label and i == 0 and j == 0 else "",
                    s=4,
                    alpha=0.2
                )

def load_data(region, regionname, season, Folder):
    Current = np.load(f"{Folder}/Current_{region}_{season}.npy")
    Next = np.load(f"{Folder}/Next_{region}_{season}.npy")
    Previous = np.load(f"{Folder}/Previous_{region}_{season}.npy")
    
    nlon = np.size(Current,axis = 1)
    ndays = np.size(Current,axis = 0)

    maxval0=np.amax(Current)
    print(f"Red value, 90th, max: {np.percentile(Current, 90):.2e} {np.amax(Current):.2e}")

    Previous = Previous/maxval0
    Next = Next/maxval0
    print(f"Purple value, 90th, max: {np.percentile(Previous, 90):.2e} {np.amax(Previous):.2e}")

    return Current, Previous, Next, nlon, ndays


def Plot(region, regionname, season, seasontime, Folder, SavingFolder):
    print(region, season)
    lev3, lev0 = plot_levels(seasontime)
    Current, Previous, Next, nlon, ndays = load_data(region, regionname, season, Folder)

    # significant scatter plot
    block_size, significance_level = 2, 0.01 # size for scatter plot
    _, significant_points1, _, _ = perform_blockwise_ttest(Current, block_size,np.mean(Current),significance_level)
    _, significant_points2, _, _ = perform_blockwise_ttest(Next, block_size, np.mean(Next),significance_level)
    _, significant_points3, _, _ = perform_blockwise_ttest(Previous, block_size,np.mean(Previous),significance_level)

    fig, ax = plt.subplots(1,1,figsize = (7,7))
  
    X = np.linspace(-int(nlon/2),int(nlon/2),nlon)
    Y = np.linspace(-int(ndays/2),int(ndays/2),ndays)

    ## Next/Previous
    contour=plt.contourf(X, Y, Previous,lev0, cmap = 'Purples',alpha=1, extend = 'max')
    plt.contourf(X, Y, Next,lev0, cmap = 'Purples', alpha=1, extend = 'max')
    cbar2 = plt.colorbar(contour, label = 'Relative Ratio')
    cbar2.set_label('Relative Ratio', fontsize=12)

    # Current Block
    plt.contourf(X, Y, Current,lev3, cmap = 'Reds' ,alpha=0.8, extend = 'max')

    cbar1 = plt.colorbar(label='LWA')
    cbar1.set_label('LWA', fontsize=12)
    
    plot_significant_points(significant_points2, block_size, nlon, ndays, label=False)
    plot_significant_points(significant_points3, block_size, nlon, ndays, label=False)
    plot_significant_points(significant_points1, block_size, nlon, ndays, label=False)

    plt.suptitle(f"Composites of LWA, {regionname} Blocks ({season})", fontsize = 14, x = 0.45, y = 0.96)
    plt.xlabel("Relative Longitudes", fontsize = 12)
    plt.ylabel("Days", fontsize = 12)
    plt.tight_layout(pad = 1.5)
    plt.savefig(f"{SavingFolder}/{region}_{season}.png",dpi=600)

    return region, season, contour

#%% main code
basepath = os.path.expanduser("~/Github")
Folder = f"{basepath}/data/ERA5/Comp_ERA"
savingpath = f'{basepath}/plots/Comp_LWA'
os.makedirs(savingpath, exist_ok=True)

regionlist = ["Pacific", "Atlantic", "BAM"]
namelist = ["Northern Pacific", "Northern Atlantic", "Southern Pacific"]

for seasontime, seasonlist in zip(["summer", "winter"], [["JJA", "JJA", "DJF"], ["DJF", "DJF", "JJA"]]):
    for region, regionname, season in zip(regionlist, namelist, seasonlist):
        region, season, contour = Plot(region, regionname, season, seasontime, Folder, savingpath)

