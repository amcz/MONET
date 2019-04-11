from . import grids, models, monet_accessor, obs, plots, sat, util

# from .monetmodels, obs, plots, util

# __all__ = ['models', 'obs', 'plots', 'sat',
#            'verification', 'util', 'monet_accessor', 'grids']
__all__ = ['models', 'obs', 'plots', 'sat', 'util', 'monet_accessor', 'grids']


def rename_latlon(ds):
    if 'latitude' in ds.coords:
        return ds.rename({'latitude': 'lat', 'longitude': 'lon'})
    elif 'Latitude' in ds.coords:
        return ds.rename({'Latitude': 'lat', 'Longitude': 'lon'})
    elif 'Lat' in ds.coords:
        return ds.rename({'Lat': 'lat', 'Lon': 'lon'})
    else:
        return ds


def rename_to_monet_latlon(ds):
    if 'lat' in ds.coords:
        return ds.rename({'lat': 'latitude', 'lon': 'longitude'})
    elif 'Latitude' in ds.coords:
        return ds.rename({'Latitude': 'latitude', 'Longitude': 'longitude'})
    elif 'Lat' in ds.coords:
        return ds.rename({'Lat': 'latitude', 'Lon': 'longitude'})
    elif 'grid_lat' in ds.coords:
        return ds.rename({'grid_lat': 'latitude', 'grid_lon': 'longitude'})
    else:
        return ds
