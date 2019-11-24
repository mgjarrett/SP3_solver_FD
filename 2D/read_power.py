import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from subprocess import call

def make_hist(datain, n):
    dataout = np.zeros([2*n]) 
    if(datain.size == n):
        for i in range(0,n):
            dataout[2*i]   = datain[i]
            dataout[2*i+1] = datain[i]
    elif(datain.size == n+1):
        for i in range(0,n):
            dataout[2*i]   = datain[i]
            dataout[2*i+1] = datain[i+1]
    return dataout
    

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

    sp3name = "./sp3/takeda_30cm_SP3_h_1mm.h5"
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

    diffname = "./sp1/takeda_30cm_SP1_h_1mm.h5"
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
    percent_error = np.zeros([nmeth+2,nx,nz])
    plot_ref = np.zeros([nx,nz])
    k_error = np.zeros([nmeth+2])
    for im in range(0,nmeth): 
        #print "method = %s " % methods[im]
        tmpsum = 0.0
        maxerror = 0.0
        for ix in range(0,nx):
            for iz in range(0,nz):
                err = mpact_pinpow[im,ix,iz]-ref_pinpow[ix,iz]
                pin_error[im,iz,ix] = err
                percent_error[im,iz,ix] = 100.0*err
                #print " %8.4f" % err,
                tmpsum = tmpsum + err*err
                if(abs(err) > maxerror):
                    maxerror = abs(err)
            #print " \n",

        rms_error[im] = math.sqrt(tmpsum / ncells)
        max_error[im] = maxerror
        k_error[im] = mpact_keff[im] - ref_keff
        
    tmpsum = 0.0
    maxerror = 0.0
    im = nmeth
    for ix in range(0,nx):
        for iz in range(0,nz):
            err = sp3_pinpow[ix,iz]-ref_pinpow[ix,iz]
            pin_error[im,iz,ix] = err
            percent_error[im,iz,ix] = 100.0*err
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
            plot_ref[ix,iz]=ref_pinpow[ix,iz]
            pin_error[im,iz,ix] = err
            percent_error[im,iz,ix] = 100.0*err
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

    #clim_lower_bound = [-0.015,-0.005,-0.005,-0.001,-0.01,-0.10]
    #clim_upper_bound = [0.015,0.005,0.005,0.001,0.01,0.10]

    # In percent
    #clim_lower_bound = [-1.5,-0.5,-0.5,-0.1,-1.0,-5.0]
    #clim_upper_bound = [1.5,0.5,0.5,0.1,1.0,5.0]

    clim_lower_bound = [-1.5,-1.5,-1.5,-1.5,-1.5,-5.0]
    clim_upper_bound = [1.5,1.5,1.5,1.5,1.5,5.0]

    Xmesh = np.linspace(0.0,15.0,num=nx+1)
    Ymesh = np.linspace(0.0,15.0,num=nz+1)

    xticks = np.linspace(0.0,15.0,num=6)
    yticks = np.linspace(0.0,15.0,num=6)

    #font = {'family' : 'normal',
    #        'weight' : 'normal',
    #        'size'   : 16}
    font = {'size' : 16}

    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    LARGE_SIZE = 22
    #plt.rc('font', size=MEDIUM_SIZE)
    #plt.rc('axes', titlesize=MEDIUM_SIZE)
    #plt.rc('axes', labelsize=MEDIUM_SIZE)
    #plt.rc('xtick', labelsize=MEDIUM_SIZE)
    #plt.rc('ytick', labelsize=MEDIUM_SIZE)
    #plt.rc('legend', fontsize=MEDIUM_SIZE)
    #plt.rc('figure', titlesize=LARGE_SIZE)

    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=SMALL_SIZE)
    plt.rc('axes', labelsize=SMALL_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=MEDIUM_SIZE)

    ### Plot reference pin power
    fig = plt.figure(0)
    #fig.set_size_inches(9, 8, forward=True)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[0,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,plot_ref[:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    #titlestr = "Power Shape Error (%s), %s" % ('%',methods[0])
    titlestr = "Reference Power Shape"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([0.5,1.5])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    #filename = "takeda_30cm_power_ref.png"
    filename = "takeda_30cm_power_ref.eps"
    dirname = "figures/reference/"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    fig.savefig(savename, format='eps', dpi=1000)
    #fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(0)
    #fig.set_size_inches(9, 8, forward=True)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[0,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,percent_error[0,:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    #titlestr = "Power Shape Error (%s), %s" % ('%',methods[0])
    titlestr = "Power Shape Error (%), ISO TL"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([clim_lower_bound[0],clim_upper_bound[0]])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_power_error_%s.png" % methods[0]
    #filename = "takeda_30cm_power_error_%s.eps" % methods[0]
    dirname = "figures/%s/" % methods[0]
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(1)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[1,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,percent_error[1,:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    #titlestr = "Power Shape Error (%s), %s" % ('%',methods[1])
    titlestr = "Power Shape Error (%), ODD TL"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([clim_lower_bound[1],clim_upper_bound[1]])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_power_error_%s.png" % methods[1]
    #filename = "takeda_30cm_power_error_%s.eps" % methods[1]
    dirname = "figures/%s/" % methods[1]
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(2)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[2,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,percent_error[2,:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    #titlestr = "Power Shape Error (%s), %s" % ('%',methods[2])
    titlestr = "Power Shape Error (%), MOM TL"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([clim_lower_bound[2],clim_upper_bound[2]])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_power_error_%s.png" % methods[2]
    #filename = "takeda_30cm_power_error_%s.eps" % methods[2]
    dirname = "figures/%s/" % methods[2]
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(3)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[3,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,percent_error[3,:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    #titlestr = "Power Shape Error (%s), %s" % ('%',methods[3])
    titlestr = "Power Shape Error (%), 1D SN"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([clim_lower_bound[3],clim_upper_bound[3]])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_power_error_%s.png" % methods[3]
    #filename = "takeda_30cm_power_error_%s.eps" % methods[3]
    dirname = "figures/%s/" % methods[3]
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(4)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[4,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,percent_error[4,:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    titlestr = "Power Shape Error (%), SP3, 1 mm"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([clim_lower_bound[4],clim_upper_bound[4]])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_power_error_SP3_1mm.png"
    #filename = "takeda_30cm_power_error_SP3_1mm.eps"
    dirname = "figures/SP3/"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(5)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)
    #plt.pcolor(Xmesh,Ymesh,pin_error[5,:,:],cmap='coolwarm')
    plt.pcolor(Xmesh,Ymesh,percent_error[5,:,:],cmap='coolwarm')
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.axis([0,15.0,0,15.0])
    titlestr = "Power Shape Error (%), SP1, 1 mm"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Axial Height [cm]')
    plt.clim([clim_lower_bound[5],clim_upper_bound[5]])
    plt.colorbar()
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_power_error_SP1_1mm.png"
    #filename = "takeda_30cm_power_error_SP1_1mm.eps"
    dirname = "figures/SP1/"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)


    ###########################################
    # comparing axial and radial power shapes #
    ###########################################

    axial_ref = np.zeros([nz])
    radial_ref = np.zeros([nx])
    axial_power = np.zeros([nmeth+2,nz])
    radial_power = np.zeros([nmeth+2,nz])
    axial_error = np.zeros([nmeth+2,nz])
    radial_error = np.zeros([nmeth+2,nz])

    ### radial power
    for ix in range(0,nx):
        tmpsum = 0.0
        for iz in range(0,nz):
            tmpsum = tmpsum + ref_pinpow[ix,iz] 
        radial_ref[ix] = tmpsum
    radnorm = float(nx)/sum(radial_ref[:])
    radial_ref[:] = radial_ref[:]*radnorm

    ### axial power
    for iz in range(0,nz):
        tmpsum = 0.0
        for ix in range(0,nx):
            tmpsum = tmpsum + ref_pinpow[ix,iz] 
        axial_ref[iz] = tmpsum
    axnorm = float(nz)/sum(axial_ref[:])
    axial_ref[:] = axial_ref[:]*axnorm

    for im in range(0,nmeth): 
        #print "method = %s " % methods[im]

        ### radial power
        for ix in range(0,nx):
            tmpsum = 0.0
            for iz in range(0,nz):
                tmpsum = tmpsum + mpact_pinpow[im,ix,iz] 
            radial_power[im,ix] = tmpsum*radnorm
            radial_error[im,ix] = radial_power[im,ix] - radial_ref[ix]

        ### axial power
        for iz in range(0,nz):
            tmpsum = 0.0
            for ix in range(0,nx):
                tmpsum = tmpsum + mpact_pinpow[im,ix,iz] 
            axial_power[im,iz] = tmpsum*axnorm
            axial_error[im,iz] = axial_power[im,iz] - axial_ref[iz]

    ### SP3 ###
    ### radial power
    for ix in range(0,nx):
        tmpsum = 0.0
        for iz in range(0,nz):
            tmpsum = tmpsum + sp3_pinpow[ix,iz] 
        radial_power[nmeth,ix] = tmpsum*radnorm
        radial_error[nmeth,ix] = radial_power[nmeth,ix] - radial_ref[ix]

    ### axial power
    for iz in range(0,nz):
        tmpsum = 0.0
        for ix in range(0,nx):
            tmpsum = tmpsum + sp3_pinpow[ix,iz] 
        axial_power[nmeth,iz] = tmpsum*axnorm
        axial_error[nmeth,iz] = axial_power[nmeth,iz] - axial_ref[iz]

    ### SP1 ###
    ### radial power
    for ix in range(0,nx):
        tmpsum = 0.0
        for iz in range(0,nz):
            tmpsum = tmpsum + diffusion_pinpow[ix,iz] 
        radial_power[nmeth+1,ix] = tmpsum*radnorm
        radial_error[nmeth+1,ix] = radial_power[nmeth+1,ix] - radial_ref[ix]

    ### axial power
    for iz in range(0,nz):
        tmpsum = 0.0
        for ix in range(0,nx):
            tmpsum = tmpsum + diffusion_pinpow[ix,iz] 
        axial_power[nmeth+1,iz] = tmpsum*axnorm
        axial_error[nmeth+1,iz] = axial_power[nmeth+1,iz] - axial_ref[iz]

    ### plot the power shapes ###
    nxhist = 30
    nzhist = 30
    plot_x = np.zeros([nxhist])
    plot_y = np.zeros([nzhist])
    plot_radial_power = np.zeros([nmeth+3,nxhist])
    plot_radial_diff  = np.zeros([nmeth+2,nxhist])
    plot_axial_power = np.zeros([nmeth+3,nzhist])
    plot_axial_diff  = np.zeros([nmeth+2,nzhist])
    plot_x = make_hist(Xmesh,nx)
    plot_y = make_hist(Ymesh,nz)

    for im in range(0,nmeth+2):
        #diff = radial_error[im,:]
        diff = 100.0*radial_error[im,:]/radial_ref[:]
        #plot_diff[im,:] = 100.0*make_hist(diff,nx)
        plot_radial_diff[im,:] = make_hist(diff,nx)
        plot_radial_power[im+1,:] = make_hist(radial_power[im,:],nx)

        #diff = axial_error[im,:]
        diff = 100.0*axial_error[im,:]/axial_ref[:]
        #plot_diff[im,:] = 100.0*make_hist(diff,nz)
        plot_axial_diff[im,:] = make_hist(diff,nz)
        plot_axial_power[im+1,:] = make_hist(axial_power[im,:],nx)

    ### Errors ###

    fig = plt.figure(6)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)

    isotlline = plt.plot(plot_x[:],plot_radial_diff[0,:], 'g-', linewidth=2.0, label='P3 ISO TL')
    isoxsline = plt.plot(plot_x[:],plot_radial_diff[1,:], 'b:', linewidth=2.0, label='P3 ANISO TL')
    sp3line = plt.plot(plot_x[:],plot_radial_diff[nmeth,:], 'r--', linewidth=2.0, label='SP3')
    #sp1line = plt.plot(plot_x[:],plot_radial_diff[nmeth+1,:], '-', color="slateblue", linewidth=2.0, label='SP3')

    titlestr = "Axially Integrated Radial Power"
    #plt.title(titlestr)
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Radial Power Error (%)')
    plt.xticks(xticks)
    #plt.legend(['P3 ISO TL','P3 ANISO TL','SP3','DIFF'],loc='best')
    #plt.legend(['P3 ISO TL','P3 ANISO TL','SP3'],loc='best')
    plt.legend(['P3 Isotropic TL','P3 Anisotropic TL','SP3'],loc='best')
    plt.axis([0,nx,-1.0,1.0])
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_radial_power_error.png"
    #filename = "takeda_30cm_radial_power_error.eps"
    dirname = "figures/radial/"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(7)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)

    isotlline = plt.plot(plot_y[:],plot_axial_diff[0,:], 'g-', linewidth=2.0, label='P3 ISO TL')
    isoxsline = plt.plot(plot_y[:],plot_axial_diff[1,:], 'b:', linewidth=2.0, label='P3 ANISO TL')
    sp3line = plt.plot(plot_y[:],plot_axial_diff[nmeth,:], 'r--', linewidth=2.0, label='SP3')
    #sp1line = plt.plot(plot_y[:],plot_axial_diff[nmeth+1,:], '-', color="slateblue", linewidth=2.0, label='SP3')

    titlestr = "Radially Integrated Axial Power"
    #plt.title(titlestr)
    plt.xlabel('Axial Height [cm]')
    plt.ylabel('Axial Power Error (%)')
    plt.xticks(xticks)
    #plt.legend(['P3 ISO TL','P3 ANISO TL','SP3','DIFF'],loc='best')
    #plt.legend(['P3 ISO TL','P3 ANISO TL','SP3'],loc='best')
    plt.legend(['P3 Isotropic TL','P3 Anisotropic TL','SP3'],loc='best')
    plt.axis([0,nz,-1.0,1.5])
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_axial_power_error.png"
    #filename = "takeda_30cm_axial_power_error.eps"
    dirname = "figures/axial"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    ### Powers ###

    fig = plt.figure(8)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)

    isotlline = plt.plot(plot_x[:],plot_radial_power[0,:], 'g-', linewidth=2.0, label='P3 ISO TL')
    isoxsline = plt.plot(plot_x[:],plot_radial_power[1,:], 'b:', linewidth=2.0, label='P3 ANISO TL')
    sp3line = plt.plot(plot_x[:],plot_radial_power[nmeth,:], 'r--', linewidth=2.0, label='SP3')
    sp1line = plt.plot(plot_x[:],plot_radial_power[nmeth+1,:], '-', color="slateblue", linewidth=2.0, label='SP3')

    titlestr = "Axially Integrated Radial Power"
    #plt.title(titlestr)
    plt.xlabel('Radial Height [cm]')
    plt.ylabel('Radial Power Error (%)')
    plt.legend(['P3 ISO TL','P3 ANISO TL','SP3','DIFF'],loc='best')
    plt.xticks(xticks)
    #plt.legend(['P3 ISO TL','P3 ANISO TL','SP3'],loc='best')
    plt.axis([0,nx,0.8,1.3])
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_radial_power.png"
    #filename = "takeda_30cm_radial_power.eps"
    dirname = "figures/radial/"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    fig = plt.figure(9)
    #plt.gcf().subplots_adjust(right=0.15)
    #plt.gcf().subplots_adjust(bottom=0.15)
    #plt.rc('font', size = 16)
    #plt.rc('font', **font)

    isotlline = plt.plot(plot_y[:],plot_axial_power[0,:], 'g-', linewidth=2.0, label='P3 ISO TL')
    isoxsline = plt.plot(plot_y[:],plot_axial_power[1,:], 'b:', linewidth=2.0, label='P3 ANISO TL')
    sp3line = plt.plot(plot_y[:],plot_axial_power[nmeth,:], 'r--', linewidth=2.0, label='SP3')
    sp1line = plt.plot(plot_y[:],plot_axial_power[nmeth+1,:], '-', color="slateblue", linewidth=2.0, label='SP3')

    titlestr = "Radially Integrated Axial Power"
    #plt.title(titlestr)
    plt.xlabel('Axial Height [cm]')
    plt.ylabel('Axial Power Error (%)')
    plt.legend(['P3 ISO TL','P3 ANISO TL','SP3','DIFF'],loc='best')
    plt.xticks(xticks)
    #plt.legend(['P3 ISO TL','P3 ANISO TL','SP3'],loc='best')
    plt.axis([0,nz,0.8,1.3])
    fig.tight_layout()
    #plt.show()

    filename = "takeda_30cm_axial_power.png"
    #filename = "takeda_30cm_axial_power.eps"
    dirname = "figures/axial"
    call(["mkdir","-p",dirname])
    savename = "%s/%s" % (dirname,filename)
    #fig.savefig(savename, format='eps', dpi=1000)
    fig.savefig(savename, format='png')
    plt.close(fig)

    ref.close()
