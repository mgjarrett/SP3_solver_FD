import h5py
import numpy as np
import math
import matplotlib.pyplot as plt

def print_attrs(name, obj):
    print name
    for key, val in obj.attrs.iteritems():
        print "    %s: %s" % (key, val)

if __name__ == '__main__':

    methods = ["isotl","oddtl","momtl","sntl"]

    problem = "x30z30"

    ref_method = "moc" 

    ########################
    # Get reference powers #
    ########################

    refname = "./%s/x30y30_hyperfine.h5" % (ref_method)
    ref = h5py.File(refname, "r")
    statename = "STATE_0001"
    ppstr = "pin_powers"
    groupname = "/%s/%s" % (statename,ppstr)
    ref_pinpow = ref[groupname]
    kstr = "keff"
    kname = "/%s/%s" % (statename,kstr)
    ref_keff = ref[kname][()]

    ###########################
    # get MPACT 2D/1D results #
    ###########################

    nx = 15
    ny = 3
    nz = 15
    na = 1

    nmeth = 4
    mpact_pinpow = np.zeros([nmeth,nx,nz])
    mpact_keff = np.zeros([nmeth])
    for tl_meth in range(0,nmeth):
        meth = methods[tl_meth]
        filename = "./%s/%s_fine.h5" % (meth,problem) 
        f = h5py.File(filename, "r")
        pin_pow = f[groupname]

        for ix in range(0,nx): 
            for iz in range(0,nz): 
                mpact_pinpow[tl_meth,ix,iz] = pin_pow[0,ix,iz]

        kstr = "keff"
        eigenvalue = f[kname][()]
        mpact_keff[tl_meth] = eigenvalue

        f.close()

    ######################
    # Get 2D SP3 results #
    ######################

    sp3name = "./sp3/takeda_30cm_SP3_h_1.7mm.h5"
    sp3_pinpow = np.zeros([nx,nz])
    f = h5py.File(sp3name, "r")
    pin_pow = f[groupname]

    for ix in range(0,nx): 
        for iz in range(0,nz): 
            sp3_pinpow[ix,iz] = pin_pow[ix,0,iz]
    sp3_keff = f[kname][()]

    f.close()

    ######################
    # Get 2D SP1 results #
    ######################

    diffname = "./sp1/takeda_30cm_SP1_h_1.7mm.h5"
    diffusion_pinpow = np.zeros([nx,nz])
    f = h5py.File(diffname, "r")
    pin_pow = f[groupname]

    for ix in range(0,nx): 
        for iz in range(0,nz): 
            diffusion_pinpow[ix,iz] = pin_pow[ix,0,iz]
    diffusion_keff = f[kname][()]

    f.close()

    #outfile = "hom_refl_p3tl.pinpow" 
    #out = open(outfile,"w")
    #for a in range (0,na):
    #    for z in range (0,nz):
    #        #out.write("Slice z = %i, assembly = %i \n" % (z,a))
    #        for x in range (0,nx):
    #            for y in range (0,ny):
    #                powstring = "%.4f " % pinpow[x,y,z,a]
    #                out.write(powstring)
    #            out.write("\n")
    #        out.write("--- \n")
    #    #out.write("\n")
    #                #print "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f" % (dataset[:,x,z,a]) 
    #            #for y in ny:
    #out.close()

    #####################
    #  compare results  #
    #####################

    ncells = nx*nz
    rms_error = np.zeros([nmeth+2])
    max_error = np.zeros([nmeth+2])
    pin_error = np.zeros([nmeth+2,nx,nz])
    k_error = np.zeros([nmeth+2])
    for im in range(0,nmeth): 
        tmpsum = 0.0
        maxerror = 0.0
        for ix in range(0,nx):
            for iz in range(0,nz):
                err = mpact_pinpow[im,ix,iz]-ref_pinpow[ix,iz]
                pin_error[im,ix,iz] = err
                tmpsum = tmpsum + err*err
                if(abs(err) > maxerror):
                    maxerror = abs(err)

        rms_error[im] = math.sqrt(tmpsum / ncells)
        max_error[im] = maxerror
        k_error[im] = mpact_keff[im] - ref_keff
        
    tmpsum = 0.0
    maxerror = 0.0
    im = nmeth
    for ix in range(0,nx):
        for iz in range(0,nz):
            err = sp3_pinpow[ix,iz]-ref_pinpow[ix,iz]
            pin_error[im,ix,iz] = err
            tmpsum = tmpsum + err*err
            if(abs(err) > maxerror):
                maxerror = abs(err)

    rms_error[im] = math.sqrt(tmpsum / ncells)
    max_error[im] = maxerror
    k_error[im] = sp3_keff - ref_keff

    tmpsum = 0.0
    maxerror = 0.0
    im = nmeth+1
    for ix in range(0,nx):
        for iz in range(0,nz):
            err = diffusion_pinpow[ix,iz]-ref_pinpow[ix,iz]
            pin_error[im,ix,iz] = err
            tmpsum = tmpsum + err*err
            if(abs(err) > maxerror):
                maxerror = abs(err)

    rms_error[im] = math.sqrt(tmpsum / ncells)
    max_error[im] = maxerror
    k_error[im] = diffusion_keff - ref_keff

    #refsum = 0.0
    #sp3sum = 0.0
    #for ix in range(0,nx):
    #    for iz in range(0,nz):
    #        refsum = refsum + ref_pinpow[ix,iz]
    #        sp3sum = sp3sum + sp3_pinpow[ix,iz]
    #print "refsum = %8.4f, sp3sum = %8.4f" % (refsum,sp3sum)

    ###############
    # print table #
    ###############

    print "       " ,
    for im in range(0,nmeth):
        print " %8s " % methods[im] ,
    print " %8s " % ("2D SP3"),
    print " %8s " % ("2D DIFF"),
    print " \n ",

    print " keff ",
    for im in range(0,nmeth+2):
        print " %8.5f " % k_error[im],
    print " \n ",

    print " RMS  ",
    for im in range(0,nmeth+2):
        print " %8.5f " % rms_error[im],
    print " \n ",

    print " MAX  ",
    for im in range(0,nmeth+2):
        print " %8.5f " % max_error[im],
    print " \n ",

    ###############     
    # plot errors #
    ###############     

    #Xmesh = np.linspace(0,15.0,num=nx+1)
    #Ymesh = np.linspace(0,15.0,num=nz+1)
    #plt.pcolor(Xmesh,Ymesh,pin_error[0,:,:],colormap=ooolwarm)

    ref.close()
