
; This ncl script calculates the pressure layer thicknesses of a hybrid coordinate system
; https://www.ncl.ucar.edu/Document/Functions/Built-in/dpres_hybrid_ccm.shtml

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
    f0 = addfile(tmp_fname, "r")


    ; ====== read pressure coefficients,pressure
    P0   = f0->P0
    hyai = f0->hyai
    hybi = f0->hybi
    ps0 = f0->PS
    

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
