from util import *
import xarray as xr

# =======================================================================

def vertical_interp(file_in, var_name, level_option, output_levels,
                    interp=2, extrap=True, P0=100000):
    '''
    Interpolates data defined at model levels to pressure (p) or
    height (z) positions

    Parameters
    ----------
    file_in : string
        The file containing the data to interpolate.
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
    return file_out.variables[var_name].get_value()


# --------------------------------------------------------------------


def isentrope_interp(file_in, var_name, output_levels, pdim):
    '''
    Interpolates data defined at model levels to isosurfaces of potential temperature

    Parameters
    ----------
    file_in : string
        The file containing the date to interpolate.
    var_name : string
        Name of variable to interpolate.
    output_levels : 1D float array
        Monotonic array giving level values at which to interpolate the variable, in Kelvin
    pdim : int
        index of the vertical level dimension of the input data

    Returns
    -------
    The variable interpolated to the given vertical positions.
    '''
    
    # call NCL script
    ncl_parent = 'isentrope_interp.ncl'
    tmpout = tmpfile()
    inputs = {'var_name':var_name, 'output_levels':output_levels, 'pdim':pdim,
              'file_in':file_in, 'file_out':tmpout}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out.rename({'LEV':'theta'})[var_name]


# --------------------------------------------------------------------


def uv2sfvpF(u, v):
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
    tmpin = tmpfile()
    uu = xr.DataArray(data=u, name='U')
    vv = xr.DataArray(data=v, name='V')
    vset = xr.merge([uu, vv])
    vset.to_netcdf(tmpin, format='NETCDF4')
   
    ncl_parent = 'uv2sfvpF.ncl'
    tmpout = tmpfile()
    inputs = {'tmp_fname':tmpin, 'file_out':tmpout}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out['SF']


# --------------------------------------------------------------------


def trop_wmo(p, t):
    '''
    Determines the level of the thermal tropopause. Simply wraps the NCL function trop_wmo:
    https://www.ncl.ucar.edu/Document/Functions/Built-in/trop_wmo.shtml

    Parameters
    ----------
    p : float array
        An array of any dimensionality containing input pressure levels. The pressure values 
        must be monotonically increasing (top-to-bottom). If multi-dimensional, it must be the 
        same size and shape as t. The level dimension must be in the rightmost position. Expected
        unit for this wrapper is hectopascals.
    t : float array
        An array of any dimensionality containing the temperatures. Missing data are not allowed. 
        The units must be degrees Kelvin. If multi-dimensional, the level dimension must be in the 
        rightmost position.

    Returns
    -------
    The returned array will be dimensioned 2 x dimsizes(u), where the 0-th element 
    of the leftmost dimension contains the stream function and the 1-th element of 
    the leftmost dimension contains the velocity potential (both in ascending 
    latitude order).
    '''
    tmpin = tmpfile()
    if(type(p) is not xr.core.dataarray.DataArray):
        p = xr.DataArray(data=p, name='p')
    if(type(t) is not xr.core.dataarray.DataArray):
        t = xr.DataArray(data=t, name='T')
    vset = xr.merge([p, t])
    vset.to_netcdf(tmpin, format='NETCDF4')
   
    ncl_parent = 'trop_wmo.ncl'
    tmpout = tmpfile()
    inputs = {'tmp_fname':tmpin, 'file_out':tmpout}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out['TROP_T']


# --------------------------------------------------------------------


def wgt_vertical_n(p, t):
    '''
    Calculates a weighted vertical average and/or sum (integral):
    https://www.ncl.ucar.edu/Document/Functions/Contributed/wgt_vertical_n.shtml

    Parameters
    ----------
    x : float array
        Array to be integrated or averaged. No missing data allowed.
    dp : float array
        Pressure thicknesses computed by dpres_hybrid_ccm or dpres_plevel. These 
        are the 'weights' associated with each level. This should have the same 
        dimensionality as x.

    Returns
    -------
    An array of of one less rank than the input x
    '''
    tmpin = tmpfile()
    if(type(p) is not xr.core.dataarray.DataArray):
        x = xr.DataArray(data=x, name='X')
    if(type(t) is not xr.core.dataarray.DataArray):
        dp = xr.DataArray(data=dp, name='dp')
    vset = xr.merge([x, dp])
    vset.to_netcdf(tmpin, format='NETCDF4')
   
    ncl_parent = 'wgt_vertrical_n.ncl'
    tmpout = tmpfile()
    inputs = {'tmp_fname':tmpin, 'file_out':tmpout}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out['VERT_AVG_X']



# --------------------------------------------------------------------


def dpres_hybrid_ccm(ps, p0, hyai, hybi):
    '''
    Calculates the pressure layer thicknesses of a hybrid coordinate system.
    Simply wraps the NCL function dpres_hybrid_ccm:
    https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_hybrid_ccm.shtml

    Parameters
    ----------
    ps : float arrays
        An array of at least 2 dimensions containing surface pressure data in Pa or hPa (mb). 
        The two rightmost dimensions must be latitude and longitude.
    p0 : float
        A scalar value equal to the surface reference pressure. Must have the same units as ps.
    hyai, hybi : float array
        A one-dimensional array equal to the hybrid A, B interface coefficients

    Returns
    -------
    The returned array will be dimensioned 2 x dimsizes(u), where the 0-th element 
    of the leftmost dimension contains the stream function and the 1-th element of 
    the leftmost dimension contains the velocity potential (both in ascending 
    latitude order).
    '''
    tmpin = tmpfile()
    ps = xr.DataArray(data=ps, name='PS')
    hyai = xr.DataArray(data=hyai, name='hyai')
    hybi = xr.DataArray(data=hybi, name='hybi')
    p0 = xr.DataArray(data=p0, name='P0')
   
    vset = xr.merge([ps, hyai, hybi, p0])
    vset.to_netcdf(tmpin, format='NETCDF4')
   
    ncl_parent = 'dpres_hybrid_ccm.ncl'
    tmpout = tmpfile()
    inputs = {'tmp_fname':tmpin, 'file_out':tmpout}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out['dp']


# --------------------------------------------------------------------


def omega_to_w(omega, p, t):
    '''
    Convert omega vertical velocity (Pa/s) to (m/s).
    Simply wraps the NCL function omega_to_w:
    https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml

    Parameters
    ----------
    omega : float arrays
        A variable of any dimensionality containing omega (Pa/s) Array size and shape must 
        match p and t
    p : float
        A variable containing pressure values (Pa). Array size and shape must match omega and t
    t : float array
        A variable containing temperature (K). Array size and shape must match omega and p

    Returns
    -------
    A double array is returned if omega is double, otherwise a float array of the same size and 
    shape as omega is returned.
    '''
    tmpin = tmpfile()
    omega = xr.DataArray(data=omega, name='OMEGA')
    p = xr.DataArray(data=p, name='P')
    t = xr.DataArray(data=t, name='T')
   
    vset = xr.merge([omega, p, t])
    vset.to_netcdf(tmpin, format='NETCDF4')
   
    ncl_parent = 'omega_to_w.ncl'
    tmpout = tmpfile()
    inputs = {'tmp_fname':tmpin, 'file_out':tmpout}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out['w']


# --------------------------------------------------------------------

def jw06_l2norm(file_in, norm_type=15):
    '''
    Computes the l2 norm (eq.14 or 15) from Jablonowski+Williamson 2006  (JW06). 
    Intended to be run on output from the JW06 baroclinic wave dycore test case
    via CESM in the "steady state" configuration

    Parameters
    ----------
    file_in : string
        The file containing the JW06 run output
    norm_type : int
        Either 14, in which case eq.14 from JW06 is computed and returned, or
        15, in which case eq.15 from JW06 is compiuted and returned

    Returns
    -------
    The l2 norm as a time series, matching the time samples of the original data
    '''
    
    # call NCL script
    ncl_parent = 'jablonowski_williamson_bw2006_l2_norm.ncl'
    tmpout = tmpfile()
    inputs = {'file_in':file_in, 'file_out':tmpout, 'norm_type':norm_type}
    file_out = call_ncl(ncl_parent, inputs, retrieve=tmpout)
    return file_out.variables['l2_norm'].values



# --------------------------------------------------------------------
