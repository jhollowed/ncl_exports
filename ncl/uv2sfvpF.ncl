
; This ncl script Computes the stream function and velocity potential via spherical harmonics given u and v on a fixed grid
; https://www.ncl.ucar.edu/Document/Functions/Built-in/uv2sfvpF-1.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
     
    if (.not. isvar("tmp_fname")) then
        print("u and v should be passed via a temporary netcdf file; check Python wrapper")
        exit()
    end if
    if (.not. isvar("file_out")) then
        print("file_in must be passed via command line")
        exit()
    end if
   
 
    ; =============================================================
    ; read u,v from file
    ; =============================================================

    print("Reading u, v")
    setfileoption("nc", "FileStructure", "Advanced")
    setfileoption("nc", "Format",  "NetCDF4Classic")
    fi = addfile(tmp_fname, "r")
    u = fi->U
    v = fi->V
    

    ; =============================================================
    ; call function
    ; =============================================================

    print("calling NCL uv2sfvpF")
    sf = uv2sfvpF(u, v)


    ; =============================================================
    ; write out
    ; =============================================================
   
    print("writing out to " + file_out) 
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c") 
    filedimdef(fo,"time",-1,True)
    fo->SF = sf
    delete(fo)
    
end
