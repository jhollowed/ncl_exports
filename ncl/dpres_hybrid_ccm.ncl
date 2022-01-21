
; This ncl script Computes the stream function and velocity potential via spherical harmonics given u and v on a fixed grid
; https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2sfvpF-1.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
     
    if (.not. isvar("tmp_fname")) then
        print("ps, p0, hyai, and hybi should be passed via a temporary netcdf file; check Python wrapper")
        exit()
    end if
    if (.not. isvar("file_out")) then
        print("file_out must be passed via command line")
        exit()
    end if
   
 
    ; =============================================================
    ; read vars from file
    ; =============================================================

    print("Reading ps, p0, hyai, hybi")
    setfileoption("nc", "FileStructure", "Advanced")
    setfileoption("nc", "Format",  "NetCDF4Classic")
    fi = addfile(tmp_fname, "r")

    ; ====== assign P0 if not provided
    if (isfilevar(fi, "P0")) then
        P0   = fi->P0                 ; Pa
    else
        print("P0 not found in input file; using 100000 Pa")
        P0 = 100000.
        P0@long_name = "reference pressure"
        P0@units     = "units"
    end if

    ; ====== read pressure coefficients, surface pressure
    hyai = fi->hyai
    hybi = fi->hybi
    ps0 = fi->PS
    

    ; =============================================================
    ; call function
    ; =============================================================

    print("calling NCL dpres_hybird_ccm")
    dp = dpres_hybrid_ccm(ps0, P0, hyai, hybi)   ;Pa [kg/(m s2)]


    ; =============================================================
    ; write out
    ; =============================================================
   
    print("writing out to " + file_out) 
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c") 
    filedimdef(fo,"time",-1,True)
    fo->dp = dp
    delete(fo)
    
end
