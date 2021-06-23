import Nio
import pdb
import wrappers
import matplotlib.pyplot as plt

p = '/scratch/cjablono_root/cjablono1/hollowed/cesm2.2/cyclone_tests/clones/cesm2.2.C192.L30.RJ12__fv3_fv_sg_adj_1800__fv3_n_sponge_30/run/cesm2.2.C192.L30.RJ12__fv3_fv_sg_adj_1800__fv3_n_sponge_30.cam.h0.0001-01-01-00000.regrid.0.5x0.5.nc'
f = Nio.open_file(p, 'wr')

# this is a hack for now; there seems to be no way to natively extract the original file location
# of a NioFile object...
f.location = p

Tz = wrappers.vertical_interp(f, 'T', 'z', [100])
pdb.set_trace()
