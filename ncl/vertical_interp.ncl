
; This ncl script interpolates data defined at model levels to 
; pressure (p) or height(z) positions

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
     
    ; =============================================================
    ; constants
    ; =============================================================
    
    gravity = 9.80616     ; gravitational acceleration in m/s^2
    

    ; =============================================================
    ; these args expected from the command line
    ; =============================================================

    if (.not. isvar("var_name")) then
        print("var_name must be passed via cli as the name of the variable to interpolate")
        exit()
    end if

    if (.not. isvar("level_option")) then
        print("level_option must be passed via cli; options are 'p' or 'z'")
        exit()
    end if
    
    if (.not. isvar("output_levels")) then
        ;print("output_levels must be passed via cli as a 1-d monotonic array in mb if \
        ;       level_option=p, or a scalar in m if level_option=z")
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

    if (.not. isvar("interp")) then
        interp = 2
    end if

    if (.not. isvar("extrap")) then
        extrap = True
    end if

    if (.not. isvar("P0")) then
        P0 = 100000.      ; Pa 
        P0mb    = P0*0.01     ; hPa
    end if


    ; =============================================================
    ; read in data
    ; =============================================================

    print("reading "+var_name+" from file")
    f = addfile(file_in, "r")
    
    ntimes = dimsizes(f->time)
    nlev = dimsizes(f->lev)
    nlat = dimsizes(f->lat)
    nlon = dimsizes(f->lon)

    var = f->$var_name$
    lev = f->lev
    lat = f->lat
    lon = f->lon
    ps = f->PS
    hyam = f->hyam
    hybm = f->hybm
    
    
    ; =============================================================
    ; do interpolation
    ; =============================================================

    print("interpolating to " + level_option + " = " + output_levels)

    if (level_option.eq."p") then
        ; ---------- interpolate in the vertical to a pressure position ----------
        var_interp = vinth2p(var, hyam, hybm, output_levels, ps, interp, P0mb, 1, extrap)

    else
        ; ---------- interpolate to a height level ----------
        if (isfilevar(f, "Z3")) then
            z3 = f->Z3
            zz = z3/gravity          ; compute geopotential height
            copyatt(zz,f->Z3)
            zz@unit = "m"
        else
            phis = f->PHIS
            t    = f->T
            q    = f->Q
            hyai = f->hyai
            hybi = f->hybi
            ; zz is top-to-bottom
            zz = cz2ccm(ps, phis, (t*(1+0.61*q)), P0, hyam(::-1), hybm(::-1), \
                                                      hyai(::-1), hybi(::-1)) 
        end if
        var_interp = wrf_interp_3d_z(var, zz, output_levels)

    end if


    ; =============================================================
    ; write out
    ; =============================================================
   
    print("writing out to " + file_out) 
    setfileoption("nc", "Format",  "NetCDF4")
    system("/bin/rm -f " + file_out) ; remove if exists
    fo = addfile(file_out, "c")
    
    filedimdef(fo,"time",-1,True)
    fo->$var_name$ = var_interp
    if(level_option.eq."p") then
        fo->ps = output_levels
    else
        fo->height = output_levels
    end if
    delete(fo)

end
