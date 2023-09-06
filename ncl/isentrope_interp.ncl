
; This ncl script interpolates data defined at model levels to isentropes

begin

    ; =============================================================
    ; these args expected from the command line
    ; =============================================================

    if (.not. isvar("var_name")) then
        print("var_name must be passed via cli as the name of the variable to interpolate")
        exit()
    end if
    
    if (.not. isvar("output_levels")) then
        ;print("output_levels must be passed via cli as a 1-d monotonic array in mb if \
        ;       level_option=p, or a scalar in m if level_option=z")
        exit()
    end if
    
    if (.not. isvar("pdim")) then
        print("pdim  must be passed via cli as an integer giving the index of the lev dimension")
        exit()
    end if

    if (.not. isvar("file_in")) then
        print("file_in must be passed via command line")
        exit()
    end if

    if (.not. isvar("file_out")) then
        print("file_in must be passed via command line")
        exit()
    end if

    ; =============================================================
    ; read in data
    ; =============================================================

    print("reading "+var_name+" from file")
    f = addfile(file_in, "r")

    var   = f->$var_name$
    lev   = f->lev
    lat   = f->lat
    p0    = f->P0
    hyam = f->hyam
    hybm = f->hybm
        
    if (isfilevar(f, "PT")) then
        ; Read the potential temperature from file if exists, assume variable name "PT"
        theta = f->PT
    else
        ; If variable "PT" is not present, then construct it from the temperature field
        ; step 1: read temperature
        ; step 2: find pressure at gridpoints
        ; step 3: calculate potential temperautre at each grid point
        T  = f->T
        P  = T                                   ; transfer meta data
        ps = f->PS
        P  = pres_hybrid_ccm (ps, p0, hya, hyb)  ; calc gridpoint pressure
        P@long_name = "pressure at each hybrid level"
        P@units     = ps@units
        theta = T                                ; transfer meta data
        theta = t*(100000/P)^0.286               ; calc potential temperature
    end if
    
    ; =============================================================
    ; do interpolation
    ; =============================================================

    print("interpolating to isentropes...")
    ityp = 1     ; do interpolation on log of pressure levels
    var_interp = int2p_n_Wrap(theta, var, output_levels, ityp, pdim)

    ; =============================================================
    ; write out
    ; =============================================================
   
    print("writing out to " + file_out) 
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c")
    
    filedimdef(fo,"time",-1,True)
    fo->$var_name$ = var_interp
    fo->isentropes = output_levels
    delete(fo)

end
