import Nio
from util import *

# =======================================================================

def vertical_interp(file_in, var_name, level_option, output_levels,
                    interp=2, extrap=True, P0=100000):
    '''
    Interpolates data defined at model levels to pressure (p) or
    height (z) positions

    Parameters
    ----------
    file_in : PyNIO NioFile object
        The file containing the date to interpolate. It is expected that the location 
        of the file is stored as a string in the object's attrubute 'location'
    var_name : string
        Name of variable to interpolate. Rightmost dimensions should be nz x nx x ny.
    level_option : char
        whether to interpolate to defined pressure or height positions
        'p' for pressure, 'z' for height 
    output_levels : 1D float array
        Monotonic array giving level values at which to interpolate the variable
        if level_option == 'p', values should be in mb/hPa
        if level option == 'z', values should be in m
    interp : int, optional
        The interpolation type. Applies to the case level_option='p' only. Defaults to 2
        1 = linear
        2 = log
        3 = log log
    extrap : bool, optional
        If false, no extrapolation is done when the pressure level is outside of the data 
        range. Applies to the case level_option='p' only. Defautls to True.
    P0 : float
        Reference surface pressure in Pa

    Returns
    -------
    The variable interpolated to the given vertical positions. Dimensionality will be the same
    as the variable, except with no nx dimension
    '''
    
    # call NCL script
    ncl_parent = 'vertical_interp.ncl'
    tmpout = tmpfile()
    inputs = {'var_name':var_name, 'level_option':level_option, 
              'output_levels':output_levels, 'file_in':file_in, 
              'file_out':tmpout, 'interp':interp, 'extrap':extrap, 
              'P0':P0}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    file_out.set_option('MaskedArrayMode', 'MaskedNever')
    return file_out.variables[var_name].get_value()


# --------------------------------------------------------------------


def uv2sfvpF(u, v)
    '''
    Computes the stream function and velocity potential via spherical harmonics 
    given u and v on a fixed grid. Simply wraps the NCL function uv2sfpF:
    https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2sfvpF-1.shtml

    Parameters
    ----------
    u,v : float arrays
        Wind components (input, arrays with two or more dimensions, rightmost two 
        dimensions must be nlat x nlon)
        -- input values must be in ascending latitude order
        -- input array must be on a global grid

    Returns
    -------
    The returned array will be dimensioned 2 x dimsizes(u), where the 0-th element 
    of the leftmost dimension contains the stream function and the 1-th element of 
    the leftmost dimension contains the velocity potential (both in ascending 
    latitude order).
    '''

    ncl_parent = 'uv2sfvpF.ncl'
    tmpout = tmpfile()
    inputs = {'u:':u, 'v':v}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    file_out.set_option('MaskedArrayMode', 'MaskedNever')
    return file_out.variables['SF'].get_value()
