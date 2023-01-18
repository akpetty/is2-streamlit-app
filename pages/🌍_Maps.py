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

#@st.cache(allow_output_mutation=True) 
#@st.experimental_memo(ttl=600, persist="disk")
def make_map_bokeh(data, variable="ice_thickness", time="Jan-01 2019", dynamic=False): 
    """Make a basic map of the data 
    
    Args: 
        data (xr.Dataset): data with all variables and timesteps
        variable (str, optional): valid variable (default to "ice_thickness")
        time (datetime or string, optional): valid time selection (default to Jan-01 2019)
        dynamic (boolean, option): return plot that dynamically updates? (default to False)
    Regards: 
        arctic_map (Holomap): map 
        
    """
    sys.stdout.write("\nMaking Arctic map...")

    # Get data for just the input variable 
    to_plot = data[variable]

    # Get colorbar limits, cmap, and clabel 
    # This is a helper function in plotting_utils.py 
    clim, cmap, clabel = get_plot_settings_by_var(variable, to_plot)

    # Get data for input time 
    # The plot will be displayed even if there's no data for that month 
    # This makes it so the time slider works 
    # If no data for that month, all values are replaced with nan and an empty map is shown
    if time != "no data": 
        to_plot = to_plot.sel(time=time)
    else: 
        to_plot = to_plot.isel(time=0)
        to_plot = to_plot.where(to_plot == -9999.) # Replace all values with nan  
    
    # Set a descriptive title
    title = "time: {0}; variable: {1}".format(time, to_plot.long_name)

    # Make map with inputs 
    # https://hvplot.holoviz.org/user_guide/Geographic_Data.html?highlight=quadmesh
    arctic_map = to_plot.hvplot.quadmesh(
        y="latitude", x="longitude", 
        clim=clim, clabel=clabel, 
        projection=ccrs.NorthPolarStereo(central_longitude=-45),
        features=["coastline"], 
        cmap=cmap,
        project=True, 
        ylim=(60,90),
        frame_width=600,
        dynamic=dynamic, 
        title=title,
        rasterize=True)
    #hv.output(widget_location="bottom") # Put time slider below the map 
    sys.stdout.write("complete!")
    return arctic_map

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
#    """Main function to create the app
#
#    Args: 
#        is2_ds (xr.Dataset): preloaded input data 
#    Return: 
#        streamlit app generated 
#    """

st.set_page_config(page_title="Mapping demos", page_icon="🌍", layout = "wide")

# -------- TITLE --------
st.title("Monthly gridded winter Arctic sea ice thickness data from ICESat-2 (IS2SITMOGR4)")


# -------- SIDEBAR INFORMATION --------

#st.sidebar.title("Mapping options")
st.sidebar.markdown("Choose your date and variable of interest, and the map will automatically update.")

is2_ds = load_data_from_aws()

# -------- USER OPTIONS FOR TIME --------
# Valid time options loaded from the dataset, then converted from string to datetime object
time_options = pd.to_datetime(is2_ds.time.values).to_pydatetime()

# Make the time slider in the streamlit sidebar 
# Every month of the year is shown even though there's no summer data 
# This has to be done in an unintuitive way because of a lack of flexibility from streamlit 
# streamlit only accepts datetime objects for min_value and max_value and timedelta objects for step
chosen_time = st.sidebar.slider(
    "Time",
    min_value=time_options[0], 
    max_value=time_options[-1],
    step=timedelta(days=31), # Approximately 1 month. timedelta object has no option for month
    format="MMM YYYY" # Only month and year are shown to the user 
).strftime("%b %Y")

# Input time to feed to the mapping function 
# Will be set to "no data" if a month is chosen on the time slider that there is no data for 
# i.e June will be set to "no data"
map_input_time = chosen_time if pd.to_datetime(chosen_time) in time_options else "no data"

# -------- USER OPTIONS FOR VARIABLE --------
# A subset of the total variable options.
# I just picked variables I thought users might be interested in seeing
var_options = [
    'freeboard_int','snow_depth_int','ice_thickness_int','ice_thickness_unc','ice_type','mean_day_of_month',
    'num_segments','sea_ice_conc','freeboard', 'ice_thickness', 'snow_depth'
]

# Descriptive long names for each of the variable options 
# i.e "sea ice thickness uncertainty" (long_name) is shown as an option instead of "ice_thickness_unc" 
# I did this instead of using the shorter variable string since its more user-friendly 
var_options_long_name = [
    is2_ds[var].long_name for var in var_options
] 

# Dict of long name : variable name to make accessing easier
var_dict = dict(zip(var_options_long_name, var_options))

# Show variable  options in streamlit sidebar 
chosen_var_A = st.sidebar.selectbox("Variable A", var_options_long_name, index=2)
# Input variable sent to the mapping function is the variable name, not the long name 
# i.e "ice_thickness_unc" instead of "ice thickness uncertainty"
map_input_var_a = var_dict[chosen_var_A]
# Lastly, the longer description of the variable is shown below the selected variable 
# I thought this was nice since it gives users a better idea of the variable 
var_description = st.sidebar.markdown('Variable A description: '+is2_ds[map_input_var_a].description)

chosen_var_B = st.sidebar.selectbox("Variable B", var_options_long_name, index=0)
map_input_var_b = var_dict[chosen_var_B]
# Lastly, the longer description of the variable is shown below the selected variable 
# I thought this was nice since it gives users a better idea of the variable 
var_description = st.sidebar.markdown('Variable B description: '+is2_ds[map_input_var_b].description)

chosen_var_C = st.sidebar.selectbox("Variable C", var_options_long_name, index=1)
map_input_var_c = var_dict[chosen_var_C]
# Lastly, the longer description of the variable is shown below the selected variable 
# I thought this was nice since it gives users a better idea of the variable 
var_description = st.sidebar.markdown('Variable C description: '+is2_ds[map_input_var_c].description)


# -------- WITH THE SELECTED OPTIONS, MAKE THE ARCTIC MAP --------
arctic_mapa = make_map_bokeh(
    data=is2_ds, # All the data 
    variable=map_input_var_a, # User-selected variable
    time=map_input_time, # User-selected time
    dynamic=True # Do not dynamically update the map 
)

# Display map in app
# Render to bokeh
# Bokeh, matplotlib, Plotly 
sys.stdout.write("\nRendering...")
arctic_mapa = hv.render(arctic_mapa, backend="bokeh")
sys.stdout.write("\nPlotting...")
st.bokeh_chart(arctic_mapa, use_container_width=False)

# -------- WITH THE SELECTED OPTIONS, MAKE THE ARCTIC MAP --------
arctic_mapb = make_map_bokeh(
    data=is2_ds, # All the data 
    variable=map_input_var_b, # User-selected variable
    time=map_input_time, # User-selected time
    dynamic=False # Do not dynamically update the map 
)

# Display map in app
# Render to bokeh
# Bokeh, matplotlib, Plotly 
arctic_mapb = hv.render(arctic_mapb, backend="bokeh")
st.bokeh_chart(arctic_mapb, use_container_width=False)

# -------- WITH THE SELECTED OPTIONS, MAKE THE ARCTIC MAP --------
arctic_mapc = make_map_bokeh(
    data=is2_ds, # All the data 
    variable=map_input_var_c, # User-selected variable
    time=map_input_time, # User-selected time
    dynamic=False # Do not dynamically update the map 
)

# Display map in app
# Render to bokeh
# Bokeh, matplotlib, Plotly 
arctic_mapc = hv.render(arctic_mapc, backend="bokeh")
st.bokeh_chart(arctic_mapc, use_container_width=False)


#if __name__ == '__main__':
    
#    main()