
; This ncl script Computes the stream function and velocity potential via spherical harmonics given u and v on a fixed grid
; https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2sfvpF-1.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
     
    if (.not. isvar("u")) then
        print("u must be passed via cli as an array with two or more dimensions, rightmost two dimensions must be nlat x nlon. Input values must be in ascending latitude order, and must be on a global grid.")
        exit()
    end if

    if (.not. isvar("v")) then
        print("v must be passed via cli as an array with two or more dimensions, rightmost two dimensions must be nlat x nlon. Input values must be in ascending latitude order, and must be on a global grid.")
        exit()
    end if


    ; =============================================================
    ; call function
    ; =============================================================

    sf = uv2sfvpF(u, v)


    ; =============================================================
    ; write out
    ; =============================================================
   
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c") 
    filedimdef(fo,"time",-1,True)
    fo->SF = sf

end
