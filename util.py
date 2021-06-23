import random
import string
import pathlib
import Nio
import subprocess
import pdb
import os

ncl = '{}/ncl'.format(pathlib.Path(__file__).parent.absolute())
tmp = '{}/tmp'.format(pathlib.Path(__file__).parent.absolute())

# =======================================================================

def tmpfile():
    '''
    Generates and names a temporary netCDF file

    Return
    ------
    A random string of 10 characters
    '''
    return tmp + '/' + ''.join(random.choice(string.ascii_lowercase) for i in range(10)) + '.nc'


def to_ncl_arr(x):
    '''
    Returns the array x as a string, formatted for NCL command line input

    e.g. the array [1, 2, 3] will return the string "(/1,2,3/)"

    Parameters
    ----------
    x : float array
        The input array
    
    Returns
    -------
    The string-formatted array in NCL syntax
    '''

    return '(/' + ','.join([str(k) for k in x]) + '/)'


def call_ncl(script, args, retrieve=None):
    '''
    Make a system call to an NCL script with command line arguments

    Parameters
    ----------
    script : string
        NCL script to execute
    args : dict
        Dictionary contining argument names (keys) and values
    retrieve : string, optional
        Name of netCDF file from which to retrieve NCL output. This will be loaded
        into a PyNio NioFile object and returned to the caller, at which point the 
        output file generated by the NCL script is deleted. If None, do nothing 
        after script call
    '''

    arg_names = args.keys()
    for i, (key,val) in enumerate(args.items()):
        # prepare all arrays present in the input args in NCL syntax
        if hasattr(val, "__len__"):
            args[key] = to_ncl_arr(val)
        # add double quotes to all strings
        if isinstance(val, str):
            args[key] = "\"{}\"".format(val)
        # get source locaitons for all PyNio file objects
        if isinstance(val, Nio.NioFile):
            args[key] = "\"{}\"".format(val.input_file)

    # make system call
    print('calling {}'.format(script))
    cmd_args = ' '.join(['\'{}={}\''.format(key, args[key]) for key in args]) + ' {}/{}'.format(ncl, script)
    os.system('ncl ' + cmd_args)

    # retrieve, clean tmp directory
    if retrieve is not None:
        print('retrieving {}'.format(retrieve))
        return Nio.open_file(retrieve, 'r')
