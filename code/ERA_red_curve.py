#%% Return Periods of Blocking Events ##
'''
Figure 3: Return Periods of Blocking Events in ERA5 and Red Noise Model
This script plots the return periods of blocking events from ERA5 data and red noise model, and the theoretical distirbution curve
'''
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import xarray as xr
import os
import seaborn as sns
import matplotlib.patches as mpatches
import BlockingDetectionFunctions 

def interval(arr):
    '''return intervals between True given a boolean array'''
    intervals = []
    true_indices = np.where(arr)[0]

    if len(true_indices) >= 2:
        gaps = np.diff(true_indices) - 1
        intervals.extend(gaps)
    return intervals


def load_phi(region, season, filepath):
    ds = xr.open_dataset(filepath)
    Lon1, Lon2, _, _ = BlockingDetectionFunctions.Region_ERA(f"{region} {season}", lat_filter=True)
    LON = ((Lon1+(Lon2-Lon1)%360/2)%360)
    phi = ds['temp_corr'].sel(lon = LON)  
    return phi.values

def load_data(nc_path):
    intervals = []
    ds = xr.open_dataset(nc_path)
    recurrence = ds.recurrence.values
    start_dates = ds['start_date'].values  # shape: (event,)
    start_dates.sort()
    time_values = ds['time'].values        # shape: (time,)

    time_to_index = {t: i for i, t in enumerate(time_values)}

    time_bool = np.zeros(len(time_values), dtype=bool)

    # Loop over start_dates and mark True where it matches a time index
    for sd in start_dates:
        if sd in time_to_index:
            time_bool[time_to_index[sd]] = True
    
    array = time_bool
    intervals = interval(array)
    method = str(ds.data)
    Pk = len(start_dates)/len(time_values)
    print(f"{method}: {int(len(start_dates))} events, {int(len(time_values))} days, Pk: {Pk}")
    return np.array(intervals), Pk, recurrence

def y_curve(phi, Pk, x):
    P = (1-phi**(x))*Pk
    y = P * ((1 - P) ** (x - 1))

    y /= y.sum()
    max_x = np.where(y==np.max(y))[0][0]
    return y, max_x

def plot(return_period1, return_period2, Pk1, Pk2, phi, region_name, N=90):
    return_period1 = return_period1[return_period1<=N] # red noise
    return_period2 = return_period2[return_period2<=N] # ERA
    x1 = np.arange(N)
    y1, max_x = y_curve(phi, Pk1, x1)
    x2 = np.arange(N)
    y2, max_x = y_curve(phi, Pk2, x2)
    c1, C1, c2, C2 = "grey", "grey", "skyblue", "k"

    fig = plt.figure(figsize=(7,5))
    sns.histplot(return_period1, bins=x1, label=f'Red noise model', stat='density', kde=False, color=c1, alpha=0.2, edgecolor=None)
    sns.histplot(return_period2, bins=x2, label=f'ERA5 data', stat='density', kde=False, color=c2, alpha=0.7, edgecolor='skyblue')
    h3, = plt.plot(x2, y2, color=C2, linewidth=2, alpha=0.5,
                label=f"Predicted recurrence constrained by ERA5\n"+
                fr"φ, α = {phi:.2g}, {Pk2:.2g}"
                )

    # Create custom legend handles with correct colors
    red_patch = mpatches.Patch(color=c1, alpha=0.2, label='Red noise model')
    era5_patch = mpatches.Patch(color=c2, alpha=0.7, label='ERA5 data')

    plt.legend(handles=[red_patch, era5_patch, h3])

    plt.yticks(np.arange(0, 0.081, 0.01))
    plt.xticks(np.arange(0, 91, 10))
    plt.ylabel("Probability Density")
    plt.xlabel("Return Period (day)")
    plt.title(f"Return Periods of {region_name}")
    plt.tight_layout(pad = 2)

    return fig

#%% Main code for red noise/ ERA5 analysis
basepath = os.path.expanduser("~/Github")
region_list = ["Pacific", "Atlantic", "BAM"]
name_list = ["Northern Pacific", "Northern Atlantic", "Southern Pacific"]
H_list = ["NH", "NH", "SH"]
season_list = ["DJF", "JJA"]

phi_path = f"{basepath}/data/Red_noise/red_noise_model" # in zenodo
red_path = f"{basepath}/data/Red_noise/BlockingEvents/ReturnPeriods"
ERA_path = f"{basepath}/data/ERA5/BlockingEvents/ReturnPeriods"

savingpath = f"{basepath}/plots/ERA_Hist"
os.makedirs(savingpath, exist_ok=True)

for region, H, name in zip(region_list, H_list, name_list):
    for season in season_list:
        ## load red noise model or ERA5
        return_period_ERA, Pk_ERA, recurrence_ERA = load_data(f"{ERA_path}/{region}_{season}.nc")
        return_period_red, Pk_red, recurrence_red = load_data(f"{red_path}/{region}_{season}.nc")
        phi = load_phi(region, season, f"{phi_path}/LWA_{H}_{season}.nc")

        # Plot the distribution
        title = f"{name} blocks, {season}"
        fig = plot(return_period_red,return_period_ERA, Pk_red, Pk_ERA, phi, title, N = 90)
        plt.savefig(f"{savingpath}/{region}_{season}.png", dpi = 600)
