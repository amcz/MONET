try:
    from pyresample.kd_tree import XArrayResamplerNN
    from pyresample.geometry import SwathDefinition, AreaDefinition
    has_pyresample = True
except ImportError:
    print('PyResample not installed.  Some functionality will be lost')
    has_pyresample = False


def _ensure_swathdef_compatability(defin):
    """ensures the SwathDefinition is compatible with XArrayResamplerNN.

    Parameters
    ----------
    defin : pyresample SwathDefinition
        a pyresample.geometry.SwathDefinition instance

    Returns
    -------
    type
        Description of returned object.

    """
    import xarray as xr
    if isinstance(defin.lons, xr.DataArray):
        return defin  # do nothing
    else:
        defin.lons = xr.DataArray(defin.lons, dims=['y', 'x']).chunk()
        defin.lats = xr.DataArray(defin.lons, dims=['y', 'x']).chunk()
        return defin


def _check_swath_or_area(defin):
    """Checks for a SwathDefinition or AreaDefinition. If AreaDefinition do
    nothing else ensure compatability with XArrayResamplerNN

    Parameters
    ----------
    defin : pyresample SwathDefinition or AreaDefinition
        Description of parameter `defin`.

    Returns
    -------
    pyresample.geometry
        SwathDefinition or AreaDefinition

    """
    try:
        if isinstance(defin, SwathDefinition):
            newswath = _ensure_swathdef_compatability(defin)
        elif isinstance(defin, AreaDefinition):
            newswath = defin
        else:
            raise RuntimeError
    except RuntimeError:
        print('grid definition must be a pyresample SwathDefinition or '
              'AreaDefinition')
        return
    return newswath


def _reformat_resampled_data(orig, new, target_grid):
    """reformats the resampled data array filling in coords, name and attrs .

    Parameters
    ----------
    orig : xarray.DataArray
        original input DataArray.
    new : xarray.DataArray
        resampled xarray.DataArray
    target_grid : pyresample.geometry
        target grid is the target SwathDefinition or AreaDefinition

    Returns
    -------
    xarray.DataArray
        reformated xarray.DataArray

    """
    target_lon, target_lat = target_grid.get_lonlats_dask()
    new.name = orig.name
    new['latitude'] = (('y', 'x'), target_lat)
    new['longitude'] = (('y', 'x'), target_lon)
    new.attrs['area'] = target_grid
    return new


def resample_stratify(da, levels, vertical, axis=1):
    import stratify
    import xarray as xr
    result = stratify.interpolate(
        levels, vertical.chunk(), da.chunk(), axis=axis)
    dims = da.dims
    out = xr.DataArray(result, dims=dims)
    for i in dims:
        if i != 'z':
            out[i] = da[i]
    out.attrs = da.attrs.copy()
    if len(da.coords) > 0:
        for i in da.coords:
            if i != 'z':
                out.coords[i] = da.coords[i]
    return out


def resample_xesmf(source_da, target_da, cleanup=False, **kwargs):
    import xesmf as xe
    regridder = xe.Regridder(source_da, target_da, **kwargs)
    if cleanup:
        regridder.clean_weight_file()
    return regridder(source_da)


def resample_nearest_neighbor_pyresample_dask(source_da,
                                              target_da,
                                              radius_of_influence=1e6):
    import dask
    import xarray as xr
    from numpy import empty_like

    # stack the dims so it is 'x','y','temporary_dim' if len(dim) > 2
    stacked_source_da, stack = _stack_dims_xarray(source_da)
    # print(stacked_source_da)
    rnnp = dask.delayed(
        resample_nearest_neighbor_pyresample(
            stacked_source_da,
            target_da,
            radius_of_influence=radius_of_influence))
    # find correct shape
    if stack:
        len_stacked = len(stacked_source_da.temporary_dim)
        target_shape = (len(target_da.y), len(target_da.x), len_stacked)
        darray = dask.array.from_delayed(
            rnnp, target_shape, dtype=source_da.dtype)
    else:
        darray = dask.array.from_delayed(
            rnnp, target_da.shape, dtype=source_da.dtype)

    # check if stacked
    if stack:
        out_array = xr.DataArray(darray, dims=('y', 'x', 'temporary_dim'))
        out_array.coords['temporary_dim'] = stacked_source_da.coords[
            'temporary_dim']
    else:
        out_array = xr.DataArray(darray, dims=('y', 'x'))

    for i in target_da.coords.keys():
        out_array.coords[i] = target_da.coords[i]

    print(out_array)
    # unstack
    if stack:
        out_array = out_array.unstack('temporary_dim')

    for i in source_da.attrs.keys():
        out_array.attrs[i] = source_da.attrs[i]
    out_array.name = source_da.name
    return out_array.transpose(*source_da.dims)


def _stack_dims_xarray(da):
    import pandas as pd
    dims = pd.Series(da.dims)
    stack_dims = dims.loc[~dims.isin(['x', 'y'])].tolist()
    if len(stack_dims) == 0:
        return da, False
    else:
        return da.stack(temporary_dim=stack_dims), True


def resample_nearest_neighbor_pyresample(source_da,
                                         target_da,
                                         radius_of_influence=1e6):
    from pyresample import geometry
    from pyresample import image
    from numpy import nan
    swath_source = geometry.SwathDefinition(
        lats=source_da.latitude.values, lons=source_da.longitude.values)
    swath_target = geometry.SwathDefinition(
        lats=target_da.latitude.values, lons=target_da.longitude.values)
    if len(source_da.dims) > 2:
        con = image.ImageContainerNearest(
            source_da.transpose('y', 'x', 'temporary_dim').data,
            swath_source,
            radius_of_influence=radius_of_influence,
            fill_value=nan)
    else:
        con = image.ImageContainerNearest(
            source_da.data,
            swath_source,
            radius_of_influence=radius_of_influence,
            fill_value=nan)
    return con.resample(swath_target).image_data


def resample_dataset(data,
                     source_grid,
                     target_grid,
                     radius_of_influence=100e3,
                     resample_cache=None,
                     return_neighbor_info=False,
                     neighbours=1,
                     epsilon=0,
                     interp='nearest'):
    # first get the source grid definition
    try:
        if source_grid is None:
            raise RuntimeError
    except RuntimeError:
        print('Must include pyresample.gemoetry in the data.attrs area_def or '
              'area')
        return

    # check for SwathDefinition or AreaDefinition
    # if swath ensure it is xarray.DataArray and not numpy for chunking
    source_grid = _check_swath_or_area(source_grid)

    # set kwargs for XArrayResamplerNN
    kwargs = dict(
        source_geo_def=source_grid,
        target_geo_def=target_grid,
        radius_of_influence=radius_of_influence,
        neighbours=neighbours,
        epsilon=epsilon)
    if interp is 'nearest':
        resampler = XArrayResamplerNN(**kwargs)
    # else:
    # resampler = XArrayResamplerBilinear(**kwargs)

    # check if resample cash is none else assume it is a dict with keys
    # [valid_input_index, valid_output_index, index_array, distance_array]
    # else generate the data
    if resample_cache is None:
        valid_input_index, valid_output_index, index_array, distance_array = resampler.get_neighbour_info(
        )
    else:
        resampler.valid_input_index = resample_cache['valid_input_index']
        resampler.valid_output_index = resample_cache['valid_output_index']
        resampler.index_array = resample_cache['index_array']
        resampler.distance_array = resample_cache['distance_array']

    # now store the resampled data temporarily in temp
    temp = resampler.get_sample_from_neighbour_info(data)

    # reformat data from temp
    out = _reformat_resampled_data(data, temp, target_grid)
    if return_neighbor_info:
        resample_cache = dict(
            valid_input_index=valid_input_index,
            valid_output_index=valid_output_index,
            index_array=index_array,
            distance_array=distance_array)
        return out, resample_cache
    else:
        return out
