import xarray as xr
import sys

def flip_lon(ds):
    '''
    Adjust lon values to make sure they are within (-180, 180)
    Allows for easier plotting and consistency between datasets 
    (e.g. comparisons with observations)

    Args:
    ds = xarray dataset
    '''
    lon = 'lon'
    ds['_longitude_adjusted'] = xr.where(
        ds[lon] > 180,
        ds[lon] - 360,
        ds[lon])
    
    # reassign the new coords to as the main lon coords
    # and sort DataArray using new coordinate values
    ds = (
        ds
        .swap_dims({lon: '_longitude_adjusted'})
        .sel(**{'_longitude_adjusted': sorted(ds._longitude_adjusted)})
        .drop(lon))
    
    ds = ds.rename({'_longitude_adjusted': lon})

    return ds

# r0.1mm
def r01mm(ds, var, num_ensemble, units=None):
    '''
    Annual days > 0.1mm pr

    Args:
    ds = xarray dataset
    var = name of variable in ds
    num_ensemble = number of ensemble members
    units = optional argument to define units (may need to be converted)
    '''
    # get longitude to -180-180 for compatibility
    ds = flip_lon(ds)

    # Handle units of input variable
    if not units:
        try:
            units = ds[var].attrs['units']
        except:
            sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n variable'units' attribute is absent from variable: {var}\n either pass argument units=<string units> or assign a units attribute to the dataset \n Valid attributes for this variable are 'mm/day' and 'kg m-2 s-1'\n Convert your data to either of these units if necessary")
    if units.lower() not in ["mm/day", "kg m-2 s-1"]:
        sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n {var} 'units' attribute: {units} is not recognized \n Valid attributes for this variable are 'mm/day' and 'kg m-2 s-1'\n Convert your data to either of these units and update the metadata")
    
    if units.lower() == "kg m-2 s-1":
        ds[var] = ds[var]*86400 # convert to mm/day
        ds[var].attrs['units'] = "mm/day"

    # Get r0.1mm
    ds_out = ds.where(ds[var]> 0.1).resample(time="1YE").count()
    ds_out = ds_out.rename({var: "r01mm"})

    # get ensemble average if > 1 member
    if num_ensemble > 1:
        if "ensemble_member" not in ds_out.dims:
             sys.exit("!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n num_ensemble > 1 but no dimension 'ensemble_member' was found \n Ensure there is a dimension ensemble_member in your ds if ensemble members > 1")
        ds_out = ds_out.mean("ensemble_member")

    # return num days r0.1mm per year
    return ds_out


# r1mm
def r1mm(ds, var, num_ensemble, units=None):
    '''
    Annual days > 1mm pr

    Args:
    ds = xarray dataset
    var = name of variable in ds
    num_ensemble = number of ensemble members
    units = optional argument to define units (may need to be converted)
    '''
    # get longitude to -180-180 for compatibility
    ds = flip_lon(ds)
    
    # Handle units of input variable
    if not units:
        try:
            units = ds[var].attrs['units']
        except:
            sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n variable'units' attribute is absent from variable: {var}\n either pass argument units=<string units> or assign a units attribute to the dataset \n Valid attributes for this variable are 'mm/day' and 'kg m-2 s-1'\n Convert your data to either of these units if necessary")
    if units.lower() not in ["mm/day", "kg m-2 s-1"]:
        sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n {var} 'units' attribute: {units} is not recognized \n Valid attributes for this variable are 'mm/day' and 'kg m-2 s-1'\n Convert your data to either of these units and update the metadata")
    
    if units.lower() == "kg m-2 s-1":
        ds[var] = ds[var]*86400 # convert to mm/day
        ds[var].attrs['units'] = "mm/day"

    # Get r1mm
    ds_out = ds.where(ds[var] > 1).resample(time="1YE").count()
    ds_out = ds_out.rename({var: "r1mm"})

    # get ensemble average if > 1 member
    if num_ensemble > 1:
        if "ensemble_member" not in ds_out.dims:
             sys.exit("!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n num_ensemble > 1 but no dimension 'ensemble_member' was found \n Ensure there is a dimension ensemble_member in your ds if ensemble members > 1")
        ds_out = ds_out.mean("ensemble_member")

    # return num days r1mm per year
    return ds_out


# r10mm
def r10mm(ds, var, num_ensemble, units=None):
    '''
    Annual days > 10mm pr
    Args:
    ds = xarray dataset
    var = name of variable in ds
    num_ensemble = number of ensemble members
    units = optional argument to define units (may need to be converted)
    '''
    # get longitude to -180-180 for compatibility
    ds = flip_lon(ds)

    # Handle units of input variable
    if not units:
        try:
            units = ds[var].attrs['units']
        except:
            sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n variable'units' attribute is absent from variable: {var}\n either pass argument units=<string units> or assign a units attribute to the dataset \n Valid attributes for this variable are 'mm/day' and 'kg m-2 s-1'\n Convert your data to either of these units if necessary")
    if units.lower() not in ["mm/day", "kg m-2 s-1"]:
        sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n {var} 'units' attribute: {units} is not recognized \n Valid attributes for this variable are 'mm/day' and 'kg m-2 s-1'\n Convert your data to either of these units and update the metadata")
    
    if units.lower() == "kg m-2 s-1":
        ds[var] = ds[var]*86400 # convert to mm/day
        ds[var].attrs['units'] = "mm/day"

    # Get r0.1mm
    ds_out = ds.where(ds[var] > 10).resample(time="1YE").count()
    ds_out = ds_out.rename({var: "r10mm"})

    # get ensemble average if > 1 member
    if num_ensemble > 1:
        if "ensemble_member" not in ds_out.dims:
             sys.exit("!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n num_ensemble > 1 but no dimension 'ensemble_member' was found \n Ensure there is a dimension ensemble_member in your ds if ensemble members > 1")
        ds_out = ds_out.mean("ensemble_member")

    # return num days r10mm per year
    return ds_out


# t40C
def t40C(ds, var, num_ensemble, units=None):
    '''
    Annual days > 40°C

    Args:
    ds = xarray dataset
    var = name of variable in ds
    num_ensemble = number of ensemble members
    units = optional argument to define units (may need to be converted)
    '''
    # get longitude to -180-180 for compatibility
    ds = flip_lon(ds)

    # Handle units of input variable
    if not units:
        try:
            units = ds[var].attrs['units']
        except:
            sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n variable'units' attribute is absent from variable: {var}\n either pass argument units=<string units> or assign a units attribute to the dataset \n Valid attributes for this variable are 'C', 'F', and 'K'\n Convert your data to either of these units if necessary")
    if units.lower() not in ["c", "f", "k"]:
        sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n {var} 'units' attribute: {units} is not recognized \n Valid attributes for this variable are 'C', 'F', and 'K'\n Convert your data to either of these units and update the metadata")
    
    if units.lower() == "k":
        ds[var] = ds[var] - 273.15 # convert to C
        ds[var].attrs['units'] = "C"
    elif units.lower() == "f":
        ds[var] = (ds[var] - 32)*(5/9) # convert to C
        ds[var].attrs['units'] = "C"


    # Get t40C
    ds_out = ds.where(ds[var] > 40).resample(time="1YE").count()
    ds_out = ds_out.rename({var: "t40C"})

    # get ensemble average if > 1 member
    if num_ensemble > 1:
        if "ensemble_member" not in ds_out.dims:
             sys.exit("!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n num_ensemble > 1 but no dimension 'ensemble_member' was found \n Ensure there is a dimension ensemble_member in your ds if ensemble members > 1")
        ds_out = ds_out.mean("ensemble_member")

    # return num days t 40°C per year
    return ds_out


# fd
def fd(ds, var, num_ensemble, units=None):
    '''
    Annual days < 0°C

    Args:
    ds = xarray dataset
    var = name of variable in ds
    num_ensemble = number of ensemble members
    units = optional argument to define units (may need to be converted)
    '''
    # get longitude to -180-180 for compatibility
    ds = flip_lon(ds)

    # Handle units of input variable
    if not units:
        try:
            units = ds[var].attrs['units']
        except:
            sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n variable'units' attribute is absent from variable: {var}\n either pass argument units=<string units> or assign a units attribute to the dataset \n Valid attributes for this variable are 'C', 'F', and 'K'\n Convert your data to either of these units if necessary")
    if units.lower() not in ["c", "f", "k"]:
        sys.exit(f"!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n {var} 'units' attribute: {units} is not recognized \n Valid attributes for this variable are 'C', 'F', and 'K'\n Convert your data to either of these units and update the metadata")
    
    if units.lower() == "k":
        ds[var] = ds[var] - 273.15 # convert to C
        ds[var].attrs['units'] = "C"
    elif units.lower() == "f":
        ds[var] = (ds[var] - 32)*(5/9) # convert to C
        ds[var].attrs['units'] = "C"

    # Get fd
    ds_out = ds.where(ds[var] < 0).resample(time="1YE").count()
    ds_out = ds_out.rename({var: "fd"})

    # get ensemble average if > 1 member
    if num_ensemble > 1:
        if "ensemble_member" not in ds_out.dims:
             sys.exit("!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n num_ensemble > 1 but no dimension 'ensemble_member' was found \n Ensure there is a dimension ensemble_member in your ds if ensemble members > 1")
        ds_out = ds_out.mean("ensemble_member")

    # return num days fd per year
    return ds_out


# simple average, any variable
def simple_average(ds, var, filenameout):
    '''
    simple time average

    Args:
    ds = xarray dataset
    var = name of variable in ds
    num_ensemble = number of ensemble members
    '''
    # get longitude to -180-180 for compatibility
    ds = flip_lon(ds)
    
    # Get simple time average
    # get ensemble average if > 1 member
    if num_ensemble > 1:
        if "ensemble_member" not in ds_out.dims:
             sys.exit("!~~~~~~~~~~~~~~~~ Error processing ~~~~~~~~~~~~~~~~! \n num_ensemble > 1 but no dimension 'ensemble_member' was found \n Ensure there is a dimension ensemble_member in your ds if ensemble members > 1")
        ds_out = ds_out.mean(["time", "ensemble_member"])
    else:
        ds_mean = ds.mean("time")
    ds_mean.to_netcdf(f"{var}_{filenameout}_A.nc")
    
    # return average ds
    return ds_mean


# Based on variable input, execute function and save one modern period file and one "A" map and one "A" preindustrial map
## Handle SSP via filename 
def get_A_maps(ds, func, var, filepath, filenameout, timepd, num_ensemble, model_id):
    '''
    Breaks datasets into periods (modern, preindustrial, future)
    Adjust periods as needed
    Be careful using this one bc there are a lot of hard coded elements

    Args:
    ds = xarray dataset
    func = one of the functions above
    var = var in ds
    filepath = path to inputs
    filenameout = output file
    timepd = "modern"/"preindustrial"/"future"
    model_id = Model # or ID (for output file naming)
    '''

    if ("ssp" in filepath) | (timepd == "future"): 
        ds = ds.sel(time=slice("2060", "2089"))
        ds_out = eval(f"{func}(ds, var, num_ensemble)") 
        ds_out.mean("time").to_netcdf(f"{filenameout}_future_{model_id}A.nc")
    elif timepd == "modern":
        modern = ds.sel(time=slice("1980", "2014"))
        modern_out = eval(f"{func}(modern, var, num_ensemble)")
        modern_out.mean("time").to_netcdf(f"{filenameout}_modern_{model_id}A.nc")
    elif timepd == "preindustrial":
        preindustrial = ds.sel(time=slice("1851", "1900"))
        pre_out = eval(f"{func}(preindustrial, var, num_ensemble)")
        pre_out.mean("time").to_netcdf(f"{filenameout}_preindustrial_{model_id}A.nc")

