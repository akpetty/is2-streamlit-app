# Import dependencies
import streamlit as st
import streamlit.components.v1 as components
import xarray as xr
import pandas as pd
import sys
from utils.read_data_utils import load_data_from_aws
from utils.plotting_utils import get_plot_settings_by_var

# Plotting dependencies 
import cartopy.crs as ccrs
import hvplot.xarray
import holoviews as hv
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

# Helps avoid some weird issues with the polar projection 
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh

# Silence warnings
import logging
logging.getLogger("param").setLevel(logging.CRITICAL)

# -------------- SET UP THE APP --------------

# Set page configuration
app_title = "ICESat-2 data dashboard"
st.set_page_config(
    page_title = app_title, # This will appear on the tab  
    page_icon = ":penguin:", 
    layout = "wide"
)
# Title the app
# I think this looks kind of ugly so I left it out 
# It's quite large 
# st.title(app_title)

def make_holomap(da,variable,title=None): 
    """ Make holomap for dashboard
    NOTE: The clabel does not work :( 
    
    Args: 
        da (xr.DataArray): data to plot 
        variable (str): variable to plot; must be valid variable name in ds 
    Returns: 
        pl (holomap): plot 
    """
    # Get colorbar limits, cmap, and clabel 
    # This is a helper function in plotting_utils.py 
    clim, cmap, clabel = get_plot_settings_by_var(variable, da)
    
    # Create plot 
    # title = "variable: "+long_name+", time: "+time_str
    pl = da.hvplot.quadmesh(y="latitude", x="longitude", 
         projection=ccrs.NorthPolarStereo(central_longitude=-45), 
         features=["coastline"], colorbar=True, clim=clim, cmap=cmap, 
         project=True, ylim=(60,90), frame_width=500, dynamic=False, title=title) 
    return(pl)

def make_map_html(data): 
    """Make a basic map of the data 
    
    Args: 
        data (xr.Dataset): data with all variables and timesteps
    Regards: 
        arctic_map (Holomap): map 
        
    """
    sys.stdout.write("\nMaking Arctic map...")

    # I need to use just a small subset of the total variable and time options 
    # This holomap takes FOREVER to load with many options 
    # I've tried tinkering with it to improve performace with no success
    var_options = [
        'freeboard_int','ice_thickness_int','ice_type'
    ]
    time_options = list(data.time.values)[:6] 

    # Make map with inputs 
    arctic_map = hv.HoloMap(
        {(i, j): make_holomap(data.sel(time=i)[j], j, title="time: "+pd.to_datetime(i).strftime("%m-%Y")+", variable: "+data[j].attrs["long_name"]) 
        for i in time_options
        for j in var_options}, 
        kdims = [hv.Dimension(('i', 'time')), hv.Dimension(('j', 'variable'))]
    ) 
    
    hv.output(widget_location="bottom") # Put time slider below the map 
    sys.stdout.write("complete!")
    return arctic_map

def main(is2_ds): 
    """Main function to create the app 

    Args: 
        is2_ds (xr.Dataset): preloaded input data 
    Return: 
        streamlit app generated 
    """

    arctic_map = make_map_html(is2_ds)

    # Display map in app
    # Render to html
    renderer = hv.renderer('bokeh')
    html = renderer.static_html(arctic_map)
    components.html(html, height=600)

if __name__ == '__main__':
    # This is run each time the app is spun up 
    # The data is cached so that it isn't reloaded any time the user changes the variable/time
    is2_ds = load_data_from_aws()
    main(is2_ds)