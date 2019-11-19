import xarray as xr
from pyresample.utils import wrap_longitudes
from scipy.io import FortranFile

try:
    import fv3grid as fg
    has_fv3grid = True
except ImportError:
    has_fv3grid = False


def open_dataset(fname, dtype='f4', res='C384', tile=1):
    """Reads the binary data for FV3-CHEM input generated by prep_chem_sources.

    Parameters
    ----------
    fname : type
        Description of parameter `fname`.
    **kwargs :
        This is the kwargs for the fv3grid. Users can define the res='C384' and
        the tile=1.
        valid values for them are:
        res = 'C768' 'C384' 'C192' 'C96' 'C48'
        tile = 1 2 3 4 5 6

    Returns
    -------
    xaray DataArray
        Description of returned object.

    """
    w = FortranFile(fname)
    a = w.read_reals(dtype=dtype)
    r = int(res[1:])
    s = a.reshape((r, r), order='F')
    if has_fv3grid:
        grid = fg.get_fv3_grid(res=res, tile=tile)
        # grid = grid.set_coords(['latitude', 'longitude', 'grid_lat', 'grid_lon'])
        grid['longitude'] = wrap_longitudes(grid.longitude)
        # grid = grid.rename({'grid_lat': 'lat_b', 'grid_lon': 'lon_b'})
        name = fname.split('.bin')[0]
        grid[name] = (('x', 'y'), s)
        return grid
    else:
        print(
            'Please install the fv3grid from https://github.com/bbakernoaa/fv3grid'
        )
        print('to gain the full capability of this dataset')
        return xr.DataArray(s, dims=('x', 'y'))


def to_prepchem_binary(data, fname='output.bin', dtype='f4'):
    """Writes to binary file for prep_chem_sources.

    Parameters
    ----------
    dset : data array
        Description of parameter `dset`.
    fname : type
        Description of parameter `fname`.
    dtype : type
        Description of parameter `dtype`.

    Returns
    -------
    type
        Description of returned object.

    """
    f = FortranFile(fname, 'w')
    f.write_record(data.astype(dtype))
    f.close()
