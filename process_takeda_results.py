## Michael Jarrett
### Process Takeda H5 output
### Generate simple region-integrated group fluxes for benchmark comparison

import h5py
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import time

if __name__ == '__main__':

    datafilename = "takeda_rodded_SP1_h_2.5mm.h5"
    f = h5py.File(datafilename, "r")
    print datafilename
    groupname_1 = "STATE_0001/flux/group1"
    groupname_2 = "STATE_0001/flux/group2"

    fluxdata_g1 = f[groupname_1]
    fluxdata_g2 = f[groupname_2]

    fuelbounds_min = [0,0,0]
    fuelbounds_max = [15,15,15]
    contbounds_min = [15,0,0]
    contbounds_max = [20,5,25]

    nx = 25
    ny = 25
    nz = 25

    ntot = nx*ny*nz

    nfuel = ( (fuelbounds_max[0]-fuelbounds_min[0])*
              (fuelbounds_max[1]-fuelbounds_min[1])*
              (fuelbounds_max[2]-fuelbounds_min[2]) )

    ncont = ( (contbounds_max[0]-contbounds_min[0])*
              (contbounds_max[1]-contbounds_min[1])*
              (contbounds_max[2]-contbounds_min[2]) )

    nrefl = ntot - nfuel - ncont

    ngrp = 2

    fuelsum = np.zeros([ngrp])
    contsum = np.zeros([ngrp])
    reflsum = np.zeros([ngrp])


    for iz in range(0,nz):
        for iy in range(0,ny):
            for ix in range(0,nx):
                # figure out where cell is
                if(ix >= fuelbounds_max[0] or iy >= fuelbounds_max[1] or iz >= fuelbounds_max[2]):
                    # is it in reflector?
                    if(ix < contbounds_min[0] or ix >= contbounds_max[0] or iy >= contbounds_max[1]):
                        reflsum[0] = reflsum[0] + fluxdata_g1[ix,iy,iz]
                        reflsum[1] = reflsum[1] + fluxdata_g2[ix,iy,iz]
                    else:
                        # it is in control
                        contsum[0] = contsum[0] + fluxdata_g1[ix,iy,iz]
                        contsum[1] = contsum[1] + fluxdata_g2[ix,iy,iz]
                else:
                    # in fuel
                    fuelsum[0] = fuelsum[0] + fluxdata_g1[ix,iy,iz]
                    fuelsum[1] = fuelsum[1] + fluxdata_g2[ix,iy,iz]

    fuelsum = fuelsum[:]/float(nfuel)
    reflsum = reflsum[:]/float(nrefl)
    contsum = contsum[:]/float(ncont)

    print "Fuel: Fast flux: %8.4e, thermal flux: %8.4e" % (fuelsum[0],fuelsum[1])
    print "Reflector: Fast flux: %8.4e, thermal flux: %8.4e" % (reflsum[0],reflsum[1])
    print "Control: Fast flux: %8.4e, thermal flux: %8.4e" % (contsum[0],contsum[1])
