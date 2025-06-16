#%% Return Periods for CESM1 Pre-industrial Control Run, RCP 8.5 Run
'''
Figure 4 (a-b) & SI: Return Periods for CESM1
This script plots the return periods of blocking events for CESM1 Pre-industrial Control Run and RCP 8.5 Run with theoretical curves
'''
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
from scipy.stats import gaussian_kde
import matplotlib.patches as mpatches
basepath = os.path.expanduser("~/Github")
Savingpath = f"{basepath}/plots/CESM_Hist/"
os.makedirs(Savingpath, exist_ok=True)

def kde (Return):
    kde = gaussian_kde(Return)
    x_vals = np.linspace(Return.min(), Return.max(), 1000)  
    kde_vals = kde(x_vals)
    max_idx = np.argmax(kde_vals)
    max_x = x_vals[max_idx]
    return max_x

def load_phi(data, region, season):
    filepath = f"{basepath}/data/CESM1/regional_lwa"
    ds = xr.open_dataset(f"{filepath}/{data}/{region}_{season}.nc")
    phi = ds['temp_corr']
    mean = ds['temp_corr']
    return phi.values, mean.values

def calculate_Pk(data, season, nevent):
    if season == 'DJF':
        season_days = 90
    elif season == 'JJA':
        season_days = 92
    if data == 'Hist': 
        nyears = 1799 
    elif data == "RCP":
        nyears = 600
    Pk = nevent/(season_days*nyears)
    print(f"{data}: season_days: {season_days}, nyears: {nyears}, nevent: {nevent}")
    return Pk, season_days

def curve_y(data, region, season, nevent):
    phi, _ = load_phi(data, region, season)
    Pk, season_days = calculate_Pk(data, season, nevent)
    x = np.arange(season_days)
    P = (1-phi**(x))*Pk
    y = P * ((1 - P) ** (x - 1))
    y /= y.sum()
    max_x = np.where(y==np.max(y))[0][0]
    print(f"{region} {season}: φ: {phi:.2g}, α: {Pk:.2g}, max: {max_x}")
    return y, max_x, season_days, phi, Pk

def Plotting_with_theo_curve(region, region_name, season):

    ## CESM1 Pre-industrial Control Run
    folder = f"{basepath}/data/CESM1/BlockingEvents/ReturnPeriods/Hist"
    ds_Hist = xr.open_dataset(f"{folder}/{region}_{season}.nc")
    y_Hist, max_x_Hist, season_days, phi_Hist, Pk_Hist = curve_y("Hist", region, season, ds_Hist.event.shape[0])
    Return_Hist = ds_Hist.return_period.values[ds_Hist.return_period.values <= season_days]
    nblock_Hist = ds_Hist.event.shape[0]*100/1799

    ## CESM1 RCP 8.5 Run
    folder = f"{basepath}/data/CESM1/BlockingEvents/ReturnPeriods/RCP"
    ds_RCP = xr.open_dataset(f"{folder}/{region}_{season}.nc")
    y_RCP, max_x_RCP, season_days, phi_RCP, Pk_RCP = curve_y("RCP", region, season, ds_RCP.event.shape[0])
    Return_RCP = ds_RCP.return_period.values[ds_RCP.return_period.values<=season_days]
    nblock_RCP = ds_RCP.event.shape[0]*100/600
    
    x = np.arange(season_days)

    ## Plot Density Function
    fig = plt.figure(figsize=[9,5])
    sns.set_style('white')
    
    # CESM1 Pre-industrial Control Run
    c1, c2 = 'darkseagreen', 'darkorange'
    BIN_Hist=np.arange(season_days)
    sns.histplot(Return_Hist, bins=BIN_Hist, stat='density', kde = False, color=c1, alpha = 0.7)
    Hist_curve, = plt.plot(x, y_Hist, color=c1, linewidth=2,
            label=fr"Predicted recurrence for CTRL run, φ, α = {phi_Hist:.2g}, {Pk_Hist:.2g}")

    # CESM1 LENS RCP8.5 Run
    BIN_RCP=np.arange(season_days)
    sns.histplot(Return_RCP, bins=BIN_RCP,  stat='density', kde = False, color = c2, alpha = 0.5)
    RCP_curve, = plt.plot(x, y_RCP, color='orange', linewidth=2,
            label=fr"Predicted recurrence for RCP 8.5 run, φ, α = {phi_RCP:.2g}, {Pk_RCP:.2g}")
    
    # Create custom legend handles with correct colors
    Hist_patch = mpatches.Patch(color=c1, alpha=0.7, label = f'CESM1 Pre-industrial Control Run (402-2200)')
    RCP_patch = mpatches.Patch(color=c2, alpha=0.5, label = f'CESM1 LENS RCP 8.5 Run (2081-2100)')

    plt.legend(handles=[Hist_patch, RCP_patch, Hist_curve, RCP_curve])

    # plt.xticks(range(0,season_days,5))
    plt.xticks(np.arange(0, 91, 10))
    plt.yticks(np.arange(0, 0.081, 0.01))
    plt.title(f'Return Period of {region_name} blocks, {season}', fontsize=14)
    plt.xlabel('Return Period (day)', fontsize=12)
    plt.ylabel('Probability Density', fontsize=12)

    plt.tight_layout(pad = 1)
    plt.savefig(f"{Savingpath}/{region}_{season}.png",dpi=600)

#%%
region_list = ["Pacific", "Atlantic", "BAM"]
region_names = ["Northern Pacific", "Northern Atlantic", "Southern Pacific"]
season_list = ["DJF", "JJA"]

for season in season_list:
    for region, region_name in zip(region_list, region_names):
        Return_CESM = Plotting_with_theo_curve(region, region_name, season)

