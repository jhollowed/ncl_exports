NCL source should be placed in /ncl

The python functions that call the NCL scripts are in wrappers.py; each new script added to /ncl needs a corresponding funciton written here (name of function should match name of associated NCL sourse file)

Functions for formatting and handling the NCL call are in util.py; these should be general and should not need modification for newly added wrappers, unless special cases come up that need to be handled
