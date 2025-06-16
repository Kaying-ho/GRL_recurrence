'''
Function code
'''
import numpy as np
import xarray as xr
from scipy.stats import gaussian_kde

def Region_ERA(region, lat_filter=True):
    if region == "Atlantic JJA":
        Lon1, Lon2, Lat1, Lat2  = 300, 360, 55, 75
    elif region == "Atlantic DJF":
        Lon1, Lon2, Lat1, Lat2  = 330, 30, 40, 60
    elif region == "Pacific JJA":
        Lon1, Lon2, Lat1, Lat2  = 150, 210, 50, 70
    elif region == "Pacific DJF":
        Lon1, Lon2, Lat1, Lat2  = 120, 180, 45, 65
    elif region == "BAM DJF":
        Lon1, Lon2, Lat1, Lat2  = 170, 230, -70, -50
    elif region == "BAM JJA":
        Lon1, Lon2, Lat1, Lat2  = 240, 300, -70, -50
    
    if lat_filter==False:
        Lat1, Lat2 = -90, 90

    return np.array([Lon1, Lon2, Lat1, Lat2])

def Region_CESM(region):
    if region == "Atlantic JJA":
        Lon1, Lon2, Lat1, Lat2  = 300, 360, 55, 75
    elif region == "Atlantic DJF":
        Lon1, Lon2, Lat1, Lat2  = 330, 30, 40, 60
    elif region == "Pacific JJA":
        Lon1, Lon2, Lat1, Lat2  = 150, 210, 55, 75
    elif region == "Pacific DJF":
        Lon1, Lon2, Lat1, Lat2  = 120, 180, 45, 65
    elif region == "BAM DJF":
        Lon1, Lon2, Lat1, Lat2  = 170, 230, -65, -45
    elif region == "BAM JJA":
        Lon1, Lon2, Lat1, Lat2  = 230, 290, -70, -50

    return np.array([Lon1, Lon2, Lat1, Lat2])

def datetime_conversion(datestamp):
    """
    Convert datestamp to datetime64
    """
    # Check the data type of datestamp and convert if needed
    print(f"Original dtype: {datestamp.dtype}")
    
    # If it's not already a datetime64 type, convert it
    if not np.issubdtype(datestamp.dtype, np.datetime64):
        print("Converting to datetime64[D]...")
        try:
            # Method 1: Try direct conversion
            start_dates_np = np.array(datestamp, dtype='datetime64[D]')
        except Exception as e:
            print(f"Direct conversion failed: {e}")
            # Method 2: Try string conversion
            try:
                start_dates_np = np.array([np.datetime64(str(d)) for d in datestamp])
            except Exception as e:
                print(f"String conversion failed: {e}")
                # Method 3: Check for datetime accessor
                if hasattr(datestamp, 'dt'):
                    start_dates_np = datestamp.dt.values
                else:
                    raise TypeError("Could not convert datestamp to datetime64")
    else:
        # If it's already datetime64, just extract the values
        print("Already datetime64 type...")
        start_dates_np = datestamp
        
    return start_dates_np

def ds_datetime_conversion(ds, date_field='time'):
    """
    Convert date field (coordinate or variable) in an xarray Dataset to datetime64
    
    Parameters:
    -----------
    ds : xarray.Dataset
        The xarray dataset containing the date field to convert
    date_field : str, default='time'
        The name of the coordinate or variable containing dates to convert
        
    Returns:
    --------
    xarray.Dataset
        The dataset with the date field converted to datetime64
    """ 
    # Make a copy to avoid modifying the original
    result_ds = ds.copy()
    
    # Check if the field exists as either a coordinate or variable
    is_coord = date_field in result_ds.coords
    is_var = date_field in result_ds.data_vars
    
    if not (is_coord or is_var):
        raise ValueError(f"'{date_field}' not found in dataset coordinates or variables")
    
    # Get the data array (works for both coords and variables)
    data_array = result_ds[date_field]
    datestamp = data_array.values
    print(f"Original dtype: {datestamp.dtype}")
    
    # If it's not already a datetime64 type, convert it
    if not np.issubdtype(datestamp.dtype, np.datetime64):
        print(f"Converting {date_field} to datetime64[ns]...")
        
        # Method 2: Try numpy conversion
        try:
            converted_dates = np.array(datestamp, dtype='datetime64[ns]')
            
            new_data_array = xr.DataArray(
                converted_dates,
                dims=data_array.dims,
                coords={k: data_array.coords[k] for k in data_array.coords if k != date_field},
                attrs=data_array.attrs
            )
            
            if is_coord:
                result_ds = result_ds.assign_coords({date_field: new_data_array})
            else:  # is_var
                result_ds[date_field] = new_data_array
                
        except Exception as e:
            print(f"Direct conversion failed: {e}")
            
            # Method 3: Try string conversion
            try:
                converted_dates = np.array([np.datetime64(str(d)) for d in datestamp])
                
                new_data_array = xr.DataArray(
                    converted_dates,
                    dims=data_array.dims,
                    coords={k: data_array.coords[k] for k in data_array.coords if k != date_field},
                    attrs=data_array.attrs
                )
                
                if is_coord:
                    result_ds = result_ds.assign_coords({date_field: new_data_array})
                else:  # is_var
                    result_ds[date_field] = new_data_array
                    
            except Exception as e:
                print(f"String conversion failed: {e}")
                raise TypeError("Could not convert date field to datetime64")
    else:
        # If it's already datetime64, no conversion needed
        print("Already datetime64 type...")
    
    return result_ds

def kde(Return):
    if len(Return) <= 1:
        return None  
    kde = gaussian_kde(Return)
    x_vals = np.linspace(Return.min(), Return.max(), 1000)  
    kde_vals = kde(x_vals)
    max_idx = np.argmax(kde_vals)
    max_x = x_vals[max_idx]
    return max_x

def stat(ds):
    mean = ds.duration.mean()
    median = ds.duration.median()
    max_return = kde(ds.return_period.values)

    ds['duration_mean'] = mean
    ds['duration_median'] = median
    ds['recurrence'] = max_return

    return ds
