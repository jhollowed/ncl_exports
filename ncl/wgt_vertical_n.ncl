
; This ncl script Calculates a weighted vertical average and/or sum (integral)..
; https://www.ncl.ucar.edu/Document/Functions/Contributed/wgt_vertical_n.shtml

begin
     
    if (.not. isvar("tmp_fname")) then
        print("x and dp should be passed via a temporary netcdf file; check Python wrapper")
        exit()
    end if
    if (.not. isvar("file_out")) then
        print("file_in must be passed via command line")
        exit()
    end if
   
 
    ; =============================================================
    ; read x,dp from file
    ; =============================================================

    print("Reading x, dp")
    setfileoption("nc", "FileStructure", "Advanced")
    setfileoption("nc", "Format",  "NetCDF4Classic")
    fi = addfile(tmp_fname, "r")
    x = fi->X
    dp = fi->dp
    

    ; =============================================================
    ; call function
    ; =============================================================

    print("calling NCL wfgt_vertical_n")
    res = trop_wmo(x, dp, 0)


    ; =============================================================
    ; write out
    ; =============================================================
   
    print("writing out to " + file_out) 
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c") 
    filedimdef(fo,"time",-1,True)
    fo->VERT_AVG_X = res
    delete(fo)
    
end
