import numpy as np
import pandas as pd

def compute_vmin_vmax(da): 
    """
    Used to find the 1st and 99th percentile of values for an input DataArray
    Makes it so the colorbar is standardized across timesteps 
    """
    vmin = np.nanpercentile(da.values, 1)
    vmax = np.nanpercentile(da.values, 99)
    return vmin, vmax 

def get_plot_settings_by_var(variable, xr_da):
    """Set colobar boundaries, colormap, and colorbar label depending on variable
    Args: 
        variable (str): valid variable name 
        xr_da (xr.DataArray): data corresponding to the input variable (not a Dataset of all variables!)
    Returns: 
        clim (tuple): colorbar limits
        cmap (str): colormap 
        clabel (str): variable long name and units 
    """
    if variable in ["ice_thickness","ice_thickness_int"]: 
        clim = (0,4)
        cmap = "viridis"
    elif variable in ["freeboard","freeboard_int"]: 
        clim = (0,0.8)
        cmap = "YlOrRd"
    elif variable in ["ice_type"]: 
        clim = (0,1)
        cmap = "YlOrRd"
    elif variable in ["snow_depth","snow_depth_int"]: 
        clim = (0,0.5)
        cmap = "inferno"
    elif variable in ["snow_density"]: 
        clim = (240,330)
        cmap = "GnBu"
    elif variable in ["sea_ice_conc"]: 
        clim = (0,1)
        cmap = "Blues_r"
    else: # Manually compute colorbar limits using vmin and vmax 
        vmin, vmax = compute_vmin_vmax(xr_da)
        clim = (vmin,vmax)
        cmap = "viridis"

    # Set colorbar label 
    if variable == "ice_type": 
        # Long name and units are really long for this variable 
        # Set them here to be shorter and look nicer in the plot 
        clabel = "Sea ice type (0 = FYI, 1 = MYI)" 
    else: 
        # Otherwise make the label say "variable longname (units)"
        clabel = xr_da.long_name + " (" + xr_da.units + ")"
    return clim, cmap, clabel

def get_winter_data(da, year_start=None, start_month="Sep", end_month="Apr", force_complete_season=False):
    """ Select data for winter seasons corresponding to the input time range 
    
    Args: 
        da (xr.Dataset or xr.DataArray): data to restrict by time; must contain "time" as a coordinate 
        year_start (str, optional): year to start time range; if you want Sep 2019 - Apr 2020, set year="2019" (default to the first year in the dataset)
        start_month (str, optional): first month in winter (default to September)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
        
    Returns: 
        da_winter (xr.Dataset or xr.DataArray): da restricted to winter seasons 
    
    """
    if year_start is None: 
        print("No start year specified. Getting winter data for first year in the dataset")
        year_start = str(pd.to_datetime(da.time.values[0]).year)
    
    start_timestep = start_month+" "+str(year_start) # mon year 
    end_timestep = end_month+" "+str(int(year_start)+1) # mon year
    winter = pd.date_range(start=start_timestep, end=end_timestep, freq="MS") # pandas date range defining winter season
    months_in_da = [mon for mon in winter if mon in da.time.values] # Just grab months if they correspond to a time coordinate in da

    if len(months_in_da) > 0: 
        if (force_complete_season == True) and (all([mon in da.time.values for mon in winter])==False): 
            da_winter = None
        else: 
            da_winter = da.sel(time=months_in_da)
    else: 
        da_winter = None
        
    return da_winter
