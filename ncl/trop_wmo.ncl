
; This ncl script Determines the level of the thermal tropopause.
; https://www.ncl.ucar.edu/Document/Functions/Built-in/trop_wmo.shtml

begin
     
    if (.not. isvar("tmp_fname")) then
        print("p (hPa) and t (K) should be passed via a temporary netcdf file; check Python wrapper")
        exit()
    end if
    if (.not. isvar("file_out")) then
        print("file_in must be passed via command line")
        exit()
    end if
   
 
    ; =============================================================
    ; read p,t from file
    ; =============================================================

    print("Reading p, t")
    setfileoption("nc", "FileStructure", "Advanced")
    setfileoption("nc", "Format",  "NetCDF4Classic")
    fi = addfile(tmp_fname, "r")
    p = fi->p
    t = fi->T
    

    ; =============================================================
    ; call function
    ; =============================================================

    print("calling NCL trop_wmo")
    res = trop_wmo(p, t, 0, False)


    ; =============================================================
    ; write out
    ; =============================================================
   
    print("writing out to " + file_out) 
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c") 
    filedimdef(fo,"time",-1,True)
    fo->TROP_T = res
    delete(fo)
    
end
