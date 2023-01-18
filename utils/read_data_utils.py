import xarray as xr
import sys
import s3fs 
import streamlit as st

# Cache the data so it isn't reloaded each time you select a different options
#@st.cache(persist=True, allow_output_mutation=True) 

# Uses st.experimental_memo to only rerun when the query changes or after 10 min.
@st.experimental_memo(ttl=600, persist="disk")
def load_data_from_aws(persist=True): 
    """
    Load ICESat-2 zarrs from s3 bucket into an xarray Dataset object 
    Simplified version of Alek's function in the IS2 Jupyter Book 
    
    Args: 
        persist (bool, optional): read the data into memory? (default to True) 
    Returns: 
        is2_ds (xr.Dataset): all the data 
    
    """
    # Input variables 
    version = 'V2'
    date_str = '201811-202204'
    local_data_path = "/data/IS2SITMOGR4/", 
    bucket_name = "icesat-2-sea-ice-us-west-2"
    
    # Path to data in S3
    s3_path = (
        's3://'+bucket_name+'/IS2SITMOGR4_'+version+ 
        '/zarr/IS2SITMOGR4_'+version+'_'+date_str+'.zarr/all/'
    )
    
    # Read data from S3 bucket
    sys.stdout.write("Loading zarr from S3 bucket...")
    s3 = s3fs.S3FileSystem(anon=True)
    store = s3fs.S3Map(
        root=s3_path, 
        s3=s3, 
        check=False
    )
    is2_ds = xr.open_zarr(
        store=store, 
        # Was getting an annoying runtime warning because the 
        # default is consolidated = True
        consolidated=False 
    )
    is2_ds = is2_ds.unify_chunks()
    sys.stdout.write("complete!")
    
    # Read data into memory 
    if persist == True:
        sys.stdout.write("\nReading data into memory...")
        is2_ds = is2_ds.persist()
        sys.stdout.write("complete!")
    return is2_ds

