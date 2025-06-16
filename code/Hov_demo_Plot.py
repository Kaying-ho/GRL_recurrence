#%% Hovmoller Diagram ##
'''
Figure 1: Hovmoller Diagram of LWA
This script generates a Hovmoller diagram of LWA for Northern Atlantic 1990 DJF with annotations for blocking events.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import xarray as xr
import cftime
import pandas as pd
import BlockingDetectionFunctions 

def date(YEAR, season, time_format):
    """
    Accommodate both cftime & datetime formats
    """
    if "cftime" in str(time_format):
        print("cftime")
        date1 = cftime.DatetimeNoLeap(YEAR, 6, 1) if season == "JJA" else cftime.DatetimeNoLeap(YEAR, 12, 1)
        date2 = cftime.DatetimeNoLeap(YEAR, 8, 31) if season == "JJA" else cftime.DatetimeNoLeap(YEAR + 1, 2, 28)
    else:
        print("datetime")
        date1 = np.datetime64(str(YEAR)+"-06-01") if season == "JJA" else np.datetime64(str(YEAR)+"-12-01")
        date2 = np.datetime64(str(YEAR)+"-08-31") if season == "JJA" else np.datetime64(str(YEAR+1)+"-02-28")
    return date1, date2

def data_filter(ds, region, season, YEAR, time_name, lat_name, lon_name, apply_lat_filter=True):
    time_format = type(ds.time.values[0])
    print(f"Time format: {time_format}")

    date1, date2  = date(YEAR, season, time_format)

    ds = ds.where((ds[time_name] >= date1) & (ds[time_name] <= date2), drop=True)
    ds = ds.sortby(time_name)

    Lon1, Lon2, Lat1, Lat2  = BlockingDetectionFunctions.Region_ERA(region+" "+season, apply_lat_filter)
    longitude, latitude = ds[lon_name], ds[lat_name]  

    region_condition = xr.where(
    ((Lon1 <= Lon2) & (Lon1 <= longitude) & (longitude <= Lon2) & (Lat1 <= latitude) & (latitude <= Lat2)) |
    ((Lon1 > Lon2) & ((Lon1 <= longitude) | (longitude <= Lon2)) & (Lat1 <= latitude) & (latitude <= Lat2)),
    True, False)

    filtered_ds = ds.where(region_condition, drop=True)

    return filtered_ds

def lon_fix(lon, lwa, event_lon, region, season):
    Lon1, Lon2, _, _ = BlockingDetectionFunctions.Region_ERA(region + " " + season)

    if lon.shape[0] != lwa.shape[1]:
        print("lon & data have mismatched shapes")
        return lon, lwa  

    if (lon.diff(dim="lon").min() < 0).item():
        print("lon is not in ascending order")

    if Lon1 > Lon2:  
        print(f"original lon: {lon}")
        shift = np.abs(lon - Lon2).argmin()+1
        lon = np.roll(lon, -shift)
        lon = (lon + 180) % 360 - 180
        lwa = np.roll(lwa, -shift, axis=1)
        print(f"lon_fix: {lon}")
        event_lon = (event_lon + 180) % 360 - 180
    return lon, lwa, event_lon

def load_ERA(basepath, H, season, region, YEAR):
    LWApath = f"/Users/kaying/Github/data/ERA5/LWA_{region}_{season}_{YEAR}.nc"
    blockpath = f"{basepath}/data/ERA5/BlockingEvents/BlockingEvents_{region}_{season}_{YEAR}.nc"

    ds_block = xr.open_dataset(blockpath)
    ds_LWA = xr.open_dataset(LWApath)

    ds_block = data_filter(ds_block, region, season, YEAR, "start_date", "event_lat", "event_lon")
    ds_LWA = data_filter(ds_LWA, region, season, YEAR, "time", "lat", "lon")
    nBlock = ds_block.event.shape[0]
    data = "ERA5 data"

    return ds_block, ds_LWA, nBlock, data

def plot_lwa_hovmoller(ds_LWA, ds_block, nBlock, YEAR, region, season, data, lon_fix, savingpath=None):
    if nBlock <= 0:
        return  # Exit early if no blocks

    print(f"{nBlock} blocks around the globe in {YEAR} {region}")

    # Zonal Mean
    ds_LWA = ds_LWA.sortby("time")
    LWA_zonal = ds_LWA.LWA.mean(dim="lat")

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [9,2]}, figsize=(7,7))
    plt.subplots_adjust(hspace=0.01)

    # Prepare data
    PlotDate = BlockingDetectionFunctions.datetime_conversion(ds_LWA.time.values)
    lon, Plot = ds_LWA.lon, LWA_zonal
    lon, Plot, event_lon = lon_fix(lon, Plot, ds_block.event_lon.values, region, season)
    levs = np.linspace(np.max(Plot)/10, np.max(Plot), 11)
    # levs = np.linspace(6e08, 6e09, 13)
    
    # Plot LWA Hovmöller
    contour = ax1.contourf(lon, PlotDate, Plot, levs, cmap='Reds')
    cbar = plt.colorbar(contour, ax=ax1, label="LWA", shrink=1)
    
    # Add blocking onset points
    BlockDate = BlockingDetectionFunctions.datetime_conversion(ds_block.start_date)
    ax1.scatter(event_lon, BlockDate, marker="x", color='black', label="Blocking Onset")

    ax1.set_xlabel("Longitude")
    ax1.yaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax1.set_yticks(PlotDate[::15])
    ax1.set_xticks(np.arange(lon[0], lon[-1]+1, 15))
    ax1.set_title(f"Hovmöller diagram of LWA, {data}")
    legend = ax1.legend(fontsize=10)
    legend.get_frame().set_alpha(0.7)

    # Prepare and plot block onset matrix
    block_array = np.zeros(len(PlotDate), dtype=bool)
    for block_date in BlockDate:

        closest_idx = np.argmin(np.abs(
        (PlotDate.values if hasattr(PlotDate, 'values') else PlotDate) - 
        (block_date.values if hasattr(block_date, 'values') else block_date)))
    
        if closest_idx < len(block_array):
            block_array[closest_idx] = True
    block_array_2d = block_array.reshape(-1, 1)

    im = ax2.imshow(
        block_array_2d, cmap='binary', aspect=0.4,
        extent=[0, 1, PlotDate[0], PlotDate[-1]], origin='lower'
    )

    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel("Blocking Onset along Time", labelpad=15)

    for date in PlotDate:
        ax2.axhline(y=date, color='gray', linewidth=0.5, alpha=0.5)

    plt.tight_layout(pad=2)

    if savingpath is not None:
        plt.savefig(f"{savingpath}/{region}_{season}.png", dpi=600)

#%% main code
basepath = os.path.expanduser("~/Github")
savingpath = f'{basepath}/plots/Hov_demo'
os.makedirs(savingpath, exist_ok=True)
YEAR, season, region, region_name, H = 1990, "DJF", "Atlantic", "Northern Atlantic", "NH"
ds_block, ds_LWA, nBlock, data = load_ERA(basepath, H, season, region, YEAR)
plot_lwa_hovmoller(ds_LWA, ds_block, nBlock, YEAR, region, season, data, lon_fix, savingpath)
