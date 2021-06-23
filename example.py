import wrappers
import Nio
import pdb

f = Nio.open_file('/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/cyclone_tests/clones/cesm2.2.C192.L30.RJ12__fv3_fv_sg_adj_1800__fv3_n_sponge_30/run/cesm2.2.C192.L30.RJ12__fv3_fv_sg_adj_1800__fv3_n_sponge_30.cam.h0.0001-01-01-00000.regrid.0.5x0.5.nc')
Tz = wrappers.vertical_interp(f, 'T', 'z', [100])
pdb.set_trace()
