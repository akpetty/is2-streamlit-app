# Import dependencies
import streamlit as st
import streamlit.components.v1 as components
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
import numpy as np
import sys
from utils.read_data_utils import load_data_from_aws
from utils.plotting_utils import get_plot_settings_by_var, get_winter_data

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


# Title the app
# I think this looks kind of ugly so I left it out 
# It's quite large 
# st.title(app_title)

#@st.cache(persist=True, allow_output_mutation=True) 
def interactive_winter_comparison_lineplot(da, years=None, title="Winter comparison", frame_width=600, frame_height=350, start_month="Sep", end_month="Apr", force_complete_season=False):
    """ Make a bokeh lineplot with markers comparing monthly mean data across winter seasons 
    
    Args: 
        da (xr.DataArray): data; must contain "time" coordinate
        years (list of str): list of years for which to plot data. 2020 would correspond to the winter season defined by start month 2020 - end month 2021 (default to all unique years in da)
        title (str, optional): title to give plot (default to "Winter comparison") 
        frame_width (int, optional): width of figure (default to 600) 
        frame_height (int, optional): height of figure (default to 350) 
        start_month (str, optional): first month in winter (default to September)
        end_month (str, optional): second month in winter; this is the following calender year after start_month (default to April)
        force_complete_season (bool, optional): require that winter season returns data if and only if all months have data? i.e. if Sep and Oct have no data, return nothing even if Nov-Apr have data? (default to False) 
        
       Returns: 
           pl (bokeh lineplot) 
        
    """
    
    if years is None: 
        years = np.unique(pd.to_datetime(da.time.values).strftime("%Y")) # Unique years in the dataset 
    
    winter_means_list = []
    for year in years:
        winter_da = get_winter_data(da, year_start=year, start_month=start_month, end_month=end_month, force_complete_season=force_complete_season) # Get data from that winter 
        if winter_da is None: # In case the user inputs a year that doesn't have data, skip this loop iteration to avoid appending None
            continue
        winter_means_list.append(winter_da)
            
    # Sort by longest --> shortest. This avoids weird issues with x axis trying to be in time order 
    winter_means_list_sorted = sorted(winter_means_list, key=lambda l: (len(l), l))[::-1]

    # Combine plots and display
    i = 0
    for da_sorted in winter_means_list_sorted: 
        winter_mean = da_sorted.mean(dim=["x","y"], keep_attrs=True) # Compute mean 
        winter_mean["time"] = pd.to_datetime(da_sorted["time"].values).strftime("%b") # Reassign the time coordinate to be just the months (Nov, Dec, ect). This allows you to easily overlay the plots on top of each other, since they share an axis
        time_str = pd.to_datetime(da_sorted.time).strftime("%Y") # Get time coordinate as string value

        pl = winter_mean.hvplot(grid=True, label=""+time_str[0]+"-"+time_str[-1], frame_width=frame_width, frame_height=frame_height, line_width=2) * winter_mean.hvplot.scatter(marker='o') # Overlay scatter plot to add markers
        if i == 0:
            pl_tot = pl
        else: 
            pl_tot *= pl 
        i+=1
        
    winters_all = pl_tot.opts(hv.opts.Layout(shared_axes=True, merge_tools=True)) # Combine lineplots into a single figure 
    winters_all.opts(title=title) # Add a title 
    winters_all.opts(legend_position='bottom_right')
    return winters_all

#def main(): 
#    """Main function to create the app #
#
#    Args: 
#        is2_ds (xr.Dataset): preloaded input data 
#    Return: 
#        streamlit app generated 
#    """

# Set page configuration
# app_title = "ICESat-2 data dashboard"
st.set_page_config(page_title="Home", page_icon=":penguin:", layout = "wide")

#st.sidebar.header("Home")

# -------- SIDEBAR INFORMATION --------
#st.sidebar.title("Home")

st.sidebar.title("Home page")
st.sidebar.markdown("Add some more ancillary information here about the data or maybe the point of this app.")

# -------- TITLE --------
st.title("Monthly gridded winter Arctic sea ice thickness data from ICESat-2 (IS2SITMOGR4)")

st.markdown("NASA's Ice, Cloud, and Land Elevation Satellite-2 ([ICESat-2](https://icesat-2.gsfc.nasa.gov/)) \
    is an advanced satellite laser altimetry system specially designed to profile Earth's fast-changing polar regions. \
    The combination of meter-scale horizontal resolution and centimeter-scale vertical precision, makes ICESat-2 extremely well \
    suited for measuring the thickness of polar sea ice - ice that forms and floats on top of the Arctic and Southern Ocean. The goal of this web app is to simply showcase our winter Arctic sea ice thickness data derived from ICESat-2 freeboards and NESOSIM snow loading since its launch in fall 2018.")


st.markdown("### Mean Inner Arctic Ocean sea ice thickness time-series")

st.markdown("Add a bit more of a description")


# This is run each time the app is spun up 
# The data is cached so that it isn't reloaded any time the user changes the variable/time
is2_ds = load_data_from_aws()

is2_ds_var=is2_ds['ice_thickness_int']

winter_lineplot = interactive_winter_comparison_lineplot(is2_ds_var, title="Arctic Ocean mean winter ICESat-2 sea ice thickness")
winter_lineplot_hv = hv.render(winter_lineplot, backend="bokeh")
st.bokeh_chart(winter_lineplot_hv, use_container_width=True)


#if __name__ == '__main__':
    
#    main()