## Michael Jarrett
### 3D SP3 solver
### For verification of 2D/1D P3 in MPACT
#########################################
# This code is a 3D diffusion and SP3 solver. The code solvers a 
# diffusion-like eigenvalue problem on a 3D Cartesian mesh. The 3D SP3 
# equations are solved in the simplified form, i.e., two coupled 
# second-order diffusion-like problems. Wielandt shift cannot be 
# applied directly to these equations, so the solution is accelerated by 
# a 3D CMFD solution, which is solved by shifted inverse power iteration.
#########################################
# The 2-group Takeda cross sections are used. Fuel, reflector, void, and 
# control rod materials are given. These cross sections are already 
# homogenized (e.g., fuel is smeared with moderator), so they are optimal 
# for this solver.
#########################################

import h5py
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
import scipy.sparse.linalg
import math
import sys
import matplotlib.pyplot as plt
import time

def calc_fiss_src(keff, phi, xsnf, xschi, nmom, igp, ng):
    src=[]
    for im in range(0,nmom):
        psi=0.0
        momsrc=0.0
        for ig in range(0,ng):
            psi+=phi[im,ig]*xsnf[ig]
        momsrc=psi*xschi[igp]/keff
        src.append(momsrc)
    return src

def calc_fiss_src_accel(phi, material, xsnf, ng, ncells, keff):
    fsrc = np.zeros([ncells])
    psi=0.0
    for ic in range(0,ncells):
        mat = material_cm[ic]
        im = ic*ngrp
        for ig in range(0,ng):
            img = im + ig
            psi+=phi[img]*xsnf[mat,ig]
        fsrc[ic] = psi/keff
    return fsrc

def calc_psi(phi, xsnf, igp, ng):
    psi=0.0
    for ig in range(0,ng):
        psi+=phi[ig]*xsnf[ig]
    return psi

def calc_inscat_src(phi, xssc, igp, gmin, gmax):
    src=0.0
    for ig in range(gmin,gmax):
        src+=phi[ig]*xssc[ig]
    src-=phi[igp]*xssc[igp]
    return src

def calc_fiss_norm(phi, xsnf, h, materials, ncells, ngrp):
    fiss_norm=0.0
    for j in range(0,ncells):
        mat = materials[j]
        for ig in range(0,ngrp):
            fiss_norm+=phi[ig,j]*xsnf[mat,ig]*h[j]
    return fiss_norm

def calc_keff(phi, phid, xsnf, oldkeff, materials, ncells):
    psid = np.zeros([ncells])
    psi  = np.zeros([ncells])
    norm = [0.0,0.0]
    for j in range(0,ncells):
        mat = materials[j]
        for ig in range(0,ngrp):
            psid[j] += phid[ig,j] * xsnf[mat,ig]
            psi[j]  += phi[ig,j]  * xsnf[mat,ig]
        norm[0]+=psid[j]*psi[j]
        norm[1]+=psi[j]*psi[j]
    keff = oldkeff * norm[1]/norm[0]
    return keff

def make_hist(x):
    nx = len(x)
    histogram = np.zeros([2*nx])
    for i in range(0,nx):
        histogram[2*i]   = x[i]
        histogram[2*i+1] = x[i]
    return histogram

def get_xs1gg(scatmat, m, ig):
    
    if(len(scatmat[mat]) > 1):
        xs1=scatmom[mat][1][ig][ig]
    else:
        xs1=0.0
    return xs1

def calc_norm_two(vec,n):
    tot=0.0
    for i in range(0,n):
        tot+=vec[i]*vec[i] 
    return math.sqrt(tot)/n

#def eval_legendre(f,x):
#    val=0.0
#    fn = np.zeros([5])
#    for g in range(0,len(f)):
#        fn[g]=f[g]
#
#    val+=fn[0]
#    val+=fn[1]*x
#    val+=fn[2]*0.5*(3.0*x*x-1.0)
#    val+=fn[3]*0.5*(5.0*x**3-3.0*x)
#    val+=fn[4]*0.125*(35.0**x**4-30.0*x**2+3.0)
#
#    return val

if __name__ == '__main__':

    #ngrp = 33 
    #nmat = 1
    #nord = 4

    ## c5g7
    #ngrp = 7
    #nmat = 8
    #nord = 1

    gmres_max_inner = 200

    SP1 = 1
    SP3 = 3
    method_order = int(sys.argv[1])

    ## Takeda
    ngrp = 2
    nmat = 4
    nord = 1

    #width_X = 60.0
    width_X = 30.0
    width_Z = 30.0
    #width_Z = 1.0
    #width_Z = 60.0
    width_Y = 1.0
    #nx = 600
    nx = int(sys.argv[2])
    #nx = 1
    ny = 1
    #nz = 1
    #nz = int(sys.argv[2])
    nz = nx
    delta_x = width_X/nx
    delta_y = width_Y/ny
    delta_z = width_Z/nz
    ncells = nx*ny*nz
    case = "takeda_2D"

    matrix_size = ncells*ngrp
    #if(matrix_size < 10):
    #    direct_inv = True
    #    lu_factor = False
    #    krylov = False
    #elif(matrix_size < 30):
    #    direct_inv = False
    #    lu_factor = True
    #    krylov = False
    #else:
    #    direct_inv = False
    #    lu_factor = False
    #    krylov = True

    linear_method = int(sys.argv[3])
    if(linear_method == 1):
        direct_inv = True
        lu_factor = False
        krylov = False
    elif(linear_method == 2):
        direct_inv = False
        lu_factor = True
        krylov = False
    elif(linear_method == 3):
        direct_inv = False
        lu_factor = False
        krylov = True

    cm_to_cart = np.zeros([ncells,3],dtype=int)
    cart_to_cm = np.zeros([nx,ny,nz],dtype=int)
    nfaces = 6
    ### face indexing
    EAST   = 0
    WEST   = 1
    NORTH  = 2
    SOUTH  = 3
    TOP    = 4
    BOTTOM = 5

    facedir = np.zeros([nfaces,3],dtype=int)
    facedir[EAST,:]   = [1,0,0]  # East
    facedir[WEST,:]   = [-1,0,0] # West
    facedir[NORTH,:]  = [0,1,0]  # North
    facedir[SOUTH,:]  = [0,-1,0] # South
    facedir[TOP,:]    = [0,0,1]  # Top
    facedir[BOTTOM,:] = [0,0,-1] # Bottom
    
    neighface = np.zeros([nfaces],dtype=int)
    neighface[EAST]   = WEST # East   -> West
    neighface[WEST]   = EAST # West   -> East
    neighface[NORTH]  = SOUTH # North  -> South
    neighface[SOUTH]  = NORTH # South  -> North
    neighface[TOP]    = BOTTOM # Top    -> Bottom
    neighface[BOTTOM] = TOP # Bottom -> Top

    neighbor = np.zeros([nfaces,ncells],dtype=int)

    h = np.zeros([nfaces,ncells])
    area = np.zeros([nfaces,ncells])
    volume = np.zeros([ncells])

    tol = 1E-6

    plotcolors = ['b-','g-','c-']

    xsdir = "Takeda"

    xst = np.zeros([nmat,ngrp,nord])
    xstr = np.zeros([nmat,ngrp,nord])
    #xscap = np.zeros([nmat,ngrp])
    xsab = np.zeros([nmat,ngrp])
    xssc = np.zeros([nmat,ngrp])
    xsf = np.zeros([nmat,ngrp])
    xsnf = np.zeros([nmat,ngrp])
    xschi = np.zeros([nmat,ngrp])
    scatmom = []
    mat_scatmom = []
    scatbounds = np.zeros([nmat,ngrp,nord,2],dtype=int)
    groups = np.zeros([ngrp],dtype=int)
    Ebound = np.zeros([ngrp],dtype=float)
    
    #phi = np.zeros([nx,ny,nz,ngrp])
    phi = np.zeros([2,ngrp,ncells])
    src = np.zeros([2,ngrp,ncells])

    material_cart = np.zeros([nx,ny,nz],dtype=int)
    material_cm = np.zeros([ncells],dtype=int)

    #plt_xst = np.zeros([nmat,ngrp])

    datadir = "data"

    #############################
    ####### READ XS FILE ########
    #############################

    #for matnum in range(0,nmat):
        #xsfile = './xs/%s/mat%i.xs' % (xsdir,matnum+1)
    xsname = "takeda-1.xsl"
    xsfile = './xs/%s/%s' % (xsdir,xsname)
    f = open(xsfile)

    lines = f.readlines()
    result = lines[1].split()
    ngrp = int(result[0])
    nmat = int(result[1])
    #nord = int(result[2])
    result = lines[2].split()

    for ig in range(0,ngrp):
        groups[ig]=ig+1
        Ebound[ig]=float(result[ig])

    u=3
    mat=0
    
    maxlines = len(lines)
    while(u < maxlines):
        result = lines[u].split()
        #print result
        if(len(result) == 0):
            u=u+1
        elif(result[0] != "XSMACRO"):
            result = lines[u].split()
        else:
            # read absorption, nu-fission, fission, xschi 
            for ig in range(0,ngrp):
                u=u+1
                result = lines[u].split()
                xsnf[mat,ig]=float(result[1])
                xsf[mat,ig]=float(result[2])
                xschi[mat,ig]=float(result[3])
                ### if first column is total absorption XS
                #xscap[mat,ig]=float(result[0])-float(result[2])
                xsab[mat,ig]=float(result[0])

            # renormalize xschi
            sumchi=sum(xschi[mat,:])
            if(sumchi > 0.0):
                xschi[mat,:]/=sumchi

            # read scattering matrix
            mat_scatmom = []
            for iord in range(0,nord):
                scatmat = []
                for ig in range(0,ngrp):
                    u=u+1
                    result = lines[u].split()
                    tmp = []
                    gmin=0
                    gmax=0
                    minset=False
                    maxset=False
                    for igp in range(0,ngrp):
                        if(not maxset):
                            tmp.append(float(result[igp])) 
                            if(minset):
                                if(float(result[igp]) == 0.0):
                                    gmax=igp
                                    maxset=True
                            else:
                                if(float(result[igp]) != 0.0):
                                    gmin=igp
                                    minset=True
                        else:
                            tmp.append(0.0)

                    if(not maxset):
                        gmax=ngrp
                    scatbounds[mat,ig,iord,0]=gmin
                    scatbounds[mat,ig,iord,1]=gmax
                    scatrow = np.asarray(tmp)
                    scatmat.append(scatrow)
                mat_scatmom.append(scatmat)
            scatmom.append(mat_scatmom)

            # calculate total cross section
            for ig in range(0,ngrp):
                scatsum=0.0
                for igp in range(0,ngrp):
                    scatsum+=scatmom[mat][0][igp][ig]
                #xst[mat,ig]=xsf[mat,ig]+xscap[mat,ig]+scatsum
                xst[mat,ig]=xsab[mat,ig]+scatsum

            mat=mat+1
        u=u+1
    f.close()
    nmat=mat

    ###########################
    # print XS that were read #
    ###########################
    outfile = "xs_output.txt"
    out = open(outfile,"w")

    for matid in range(0,nmat):
        matstr  = "Material %i \n" % matid
        headstr = " %-12s %-12s %-12s %-12s \n" % ("capture","nu-fission","fission","chi")
        out.write(matstr)
        out.write(headstr)
        # write absorption, nu-fission, fission, chi
        for ig in range(0,ngrp):
            datastr = " %12.6e %12.6e %12.6e %12.6e \n" % (xsab[matid,ig],xsnf[matid,ig],xsf[matid,ig],xschi[matid,ig])
            out.write(datastr)
        # write scattering matrix
        for ig in range(0,ngrp):
            datastr = ""
            scatrow = scatmom[matid][0][ig]
            for igp in range(0,ngrp):
                #datastr = " %s %12.6e" % (datastr,scatmom[matid][0][ig][igp])
                datastr = "%s %12.6e" % (datastr,scatrow[igp])
            datastr = "%s \n" % datastr
            out.write(datastr)
        out.write("\n") 
    out.close()

    #for ig in range(0,ngrp):
    #    plt_xst[ig]=xst[ig,0]

    ##############################
    # Perform 3D SP3 calculation #
    ##############################

    sc_order=1 # order of scattering (1 = isotropic, 2 = linear, etc.)

    ninners=5 # to converge between phi0, phi2

    legendstr = ["Scalar flux"]

    ### continuous
    legendcont = ["Scalar flux"]

    keff=1.0
     
    # set up mesh 
    nrad=nx*ny
    if(case == "takeda_2D"):
        #matfile = "materials.out"
        #f = open(matfile,"w")
        #nfuelx=nx*3/4 ### reflector on east
        #nfuelx=nx*4/6 ### reflector on east
        nfuelx=nx*1/2 ### reflector on east
        #nfuelx=nx*10/11 ### reflector on east
        #nfuelx=nx*2/6 ### reflector on west
        #nfuelx=nx*5/8 ### reflector on west
        #nfuelx=nx ### no reflector
        nreflx=nx
        nfuely=ny 
        nfuelz=nz*1/2
        #nfuelz=nz*4/6
        #nfuelz=nz*1/6
        #nfuelz=nz
        nreflz=nz

        print "fuel length = %.4f" % (nfuelx*delta_x)

        ### setup mesh dimensions and materials and
        ### conversion between cartesian and 1D index
        for iz in range (0,nz):
            for iy in range (0,ny):
                for ix in range (0,nx):
                    ic=iz*nrad+iy*nx+ix
                    cm_to_cart[ic,0]=ix
                    cm_to_cart[ic,1]=iy
                    cm_to_cart[ic,2]=iz
                    cart_to_cm[ix,iy,iz]=ic
                    
                    h[0,ic] = delta_x
                    h[1,ic] = delta_x
                    h[2,ic] = delta_y
                    h[3,ic] = delta_y
                    h[4,ic] = delta_z
                    h[5,ic] = delta_z

                    volume[ic] = h[0,ic]*h[2,ic]*h[4,ic]

                    area[0,ic] = h[2,ic]*h[4,ic]
                    area[1,ic] = h[2,ic]*h[4,ic]
                    area[2,ic] = h[0,ic]*h[4,ic]
                    area[3,ic] = h[0,ic]*h[4,ic]
                    area[4,ic] = h[0,ic]*h[2,ic]
                    area[5,ic] = h[0,ic]*h[2,ic]

                    #if(ix >= nfuelx and iz < nfuelz):
                    #if(ix < nfuelx and iz >= nfuelz):
                    #if(ix < nfuelx and ix >= nfuely and iz < nfuelz):
                    if(ix < nfuelx and iz < nfuelz):
                        material_cart[ix,iy,iz] = 0 # fuel
                        material_cm[ic] = 0
                    else:
                        material_cart[ix,iy,iz] = 1 # reflector
                        material_cm[ic] = 1

                # end ix in (0,nx)
            # end iy in (0,ny)
        # end iz in (0,nz)

        # setup up neighbor relationships
        for iz in range (0,nz):
            for iy in range (0,ny):
                for ix in range (0,nx):
                    ic=iz*nrad+iy*nx+ix

                    for iface in range(0,nfaces):
                        ixp=ix+facedir[iface,0]
                        iyp=iy+facedir[iface,1]
                        izp=iz+facedir[iface,2]
                        if( (ixp < nx and iyp < ny and izp < nz) and 
                            (ixp >= 0 and iyp >= 0 and izp >= 0) ): 

                            icp=cart_to_cm[ixp,iyp,izp]
                        else:
                            icp=-1

                        neighbor[iface,ic]=icp

                    # end iface in (0,nfaces)
                # end ix in (0,nx)
            # end iy in (0,ny)
        # end iz in (0,nz)
                    
        #f.close()

    # set initial guess

    scalar_flx = np.zeros([ngrp,ncells])
    oldflx = np.zeros([ngrp,ncells])

    for ic in range(0,ncells):
        scalar_flx[:,ic]=1.0

    # zero initial partial currents
    #J_pos = np.zeros([p_order,ncells+1,ngrp])
    #J_neg = np.zeros([p_order,ncells+1,ngrp])
    surf_cur = np.zeros([nfaces,ngrp,ncells])

    psi = np.zeros([ncells])

    # set boundary conditions

    # infinite y, reflective/vacuum x and z
    bc = ["vacuum","reflective","reflective","reflective","vacuum","reflective"]
    # all reflective 
    #bc = ["reflective","reflective","reflective","reflective","reflective","reflective"]
    # infinite y and z, reflective/vacuum x (1D)
    #bc = ["vacuum","reflective","reflective","reflective","reflective","reflective"]
    #bc = ["reflective","vacuum","reflective","reflective","reflective","reflective"]
    #bc = ["vacuum","vacuum","reflective","reflective","reflective","reflective"]
    # infinite y and x, reflective/vacuum z (1D)
    #bc = ["reflective","reflective","reflective","reflective","reflective","vacuum"]
    #bc = ["reflective","reflective","reflective","reflective","vacuum","reflective"]

    fiss_norm = calc_fiss_norm(scalar_flx[:,:],xsnf,volume,material_cm,ncells,ngrp)
    fissdiff = 1.0

    #logname = "2D_SP3_takeda.log"

    fuelsize = "%icm" % (int(2.0*delta_x*nfuelx))
    if(10.0*delta_x == float(int(10.0*delta_x))):
        meshsize = "%imm" % int(10.0*delta_x)
    else:
        mesh_delta = "%2.1f" % (delta_x*10.0)
        meshsize = "%smm" % (mesh_delta)
        #print "mesh_delta = %s" % mesh_delta
        #print "mesh size  = %s" % meshsize

    if(method_order == SP3):
        method_name = "SP3"
    elif(method_order == SP1):
        method_name = "SP1"

    logname = "takeda_%s_%s_h_%s.log" % (fuelsize,method_name,meshsize)
    
    logfile = open(logname,"w")

    ####################################
    ######## Diffusion Solution ########
    ####################################

    # setup dtil
    dtil = np.zeros([2,nfaces,ngrp,ncells])
    dhat = np.zeros([nfaces,ngrp,ncells])

    for ic in range(0,ncells):
        ix=cm_to_cart[ic,0]
        iy=cm_to_cart[ic,1]
        iz=cm_to_cart[ic,2]
        matc=material_cm[ic]

        for iface in range(0,nfaces):
            hc=h[iface,ic]

            icp = neighbor[iface,ic]
            if(icp >= 0):
                matp=material_cm[icp]
                hp=h[neighface[iface],icp]

                for ig in range(0,ngrp):
                    dtil[0,iface,ig,ic]=area[iface,ic]*(2.0/3.0)/(hp*xst[matp,ig]+hc*xst[matc,ig])
                    dtil[1,iface,ig,ic]=area[iface,ic]*(18.0/35.0)/(hp*xst[matp,ig]+hc*xst[matc,ig])
                   
            else: 
                # boundary case, do something special
                #### boundary dtil #####
                for ig in range(0,ngrp):
                    if(bc[iface] == "reflective"):
                        dtil[0,iface,ig,ic]=0.0            
                        dtil[1,iface,ig,ic]=0.0            
                    elif(bc[iface] == "vacuum"):
                        #xs1p=get_xs1gg(scatmom,matc,ig)
                        #xstrp=xst[matc,ig]-xs1p
                        dtil[0,iface,ig,ic]=area[iface,ic]*0.5/(1.0+3.0*0.5*xst[matc,ig]*hc)
                        dtil[1,iface,ig,ic]=area[iface,ic]*9.0/(18.0+35.0*xst[matc,ig]*hc)

                #end ig in ngrp
            # end if boundary face
        # end iface
    # end ic

    ##### Set up zeroth order problem
    ##### Full Multigroup matrix using node-major indexing

    ###### Matrices for CMFD eigenvalue problem
    ### store as sparse

    print "Matrix size = %i " % matrix_size
    if(krylov):
        M = sparse.lil_matrix((ncells*ngrp,ncells*ngrp),dtype=float)
        F = sparse.lil_matrix((ncells*ngrp,ncells*ngrp),dtype=float)
        #M = np.zeros([ncells*ngrp,ncells*ngrp])
        #F = np.zeros([ncells*ngrp,ncells*ngrp])
    else:
        M = np.zeros([ncells*ngrp,ncells*ngrp])
        F = np.zeros([ncells*ngrp,ncells*ngrp])
        Minv = np.zeros([ncells*ngrp,ncells*ngrp])
        invMF = np.zeros([ncells*ngrp,ncells*ngrp])

    cmflx = np.ones([ncells*ngrp])
    y = np.zeros([ncells*ngrp])
    x = np.zeros([ncells*ngrp])
    x2 = np.zeros([ncells*ngrp])

    ###### Matrices for inner fixed-source problems
    if(krylov):
        A0 = sparse.lil_matrix((ncells*ngrp,ncells*ngrp),dtype=float)
        A2 = sparse.lil_matrix((ncells*ngrp,ncells*ngrp),dtype=float)
        #A0 = np.zeros([ncells*ngrp,ncells*ngrp])
        #A2 = np.zeros([ncells*ngrp,ncells*ngrp])
    else:
        A0 = np.zeros([ncells*ngrp,ncells*ngrp])
        A0inv = np.zeros([ncells*ngrp,ncells*ngrp])
        if(method_order == SP3):
            A2 = np.zeros([ncells*ngrp,ncells*ngrp])
            A2inv = np.zeros([ncells*ngrp,ncells*ngrp])

    Q0 = np.zeros([ncells*ngrp])
    Q2 = np.zeros([ncells*ngrp])

    ##### Setup coefficient matrices for fixed-source problems.        #####
    ##### These coefficients remain constant throughout the iteration. #####

    for ic in range(0,ncells):
        im = ic*ngrp
        mat = material_cm[ic]

        ### Set up transport/scattering matrix
        for ig in range(0,ngrp):
            img=im+ig
            aii = 0.0
            bii = 0.0

            #xsR0=xst[mat,ig]  ### self-scattering is subtracted later (not for SP3)
            xsR0=xst[mat,ig]-scatmom[mat][0][ig][ig]
            xsR2=xst[mat,ig]+0.8*(xst[mat,ig]-scatmom[mat][0][ig][ig])
            
            #### loop through all neighbor surfaces, add diffusion term
            for iface in range(0,nfaces):
                icp = neighbor[iface,ic]

                if(icp >= 0):
                    matp = material_cm[icp]
                    hp = h[neighface[iface],icp]
                    ipg = icp*ngrp + ig

                    A0[img,ipg] = -dtil[0,iface,ig,ic]
                    A2[img,ipg] = -dtil[1,iface,ig,ic]
                    aii = aii + dtil[0,iface,ig,ic]
                    bii = bii + dtil[1,iface,ig,ic]
                else:
                    ### boundary case, do something special
                    aii = aii + dtil[0,iface,ig,ic]
                    bii = bii + dtil[1,iface,ig,ic]

            #### add removal term
            aii = aii + xsR0*volume[ic] 
            bii = bii + xsR2*volume[ic] 
            A0[img,img]=aii
            A2[img,img]=bii

            ### set up off-diagonal scattering (not when using SP3)
            #gmin=scatbounds[mat][ig][0][0]
            #gmax=scatbounds[mat][ig][0][1]
            #for igp in range(gmin,gmax):
            #    A0[img,im+igp]=A0[img,im+igp]-scatmom[mat][0][ig][igp]*volume[ic]

            ### setup fission matrix for CMFD (also doesn't change during iteration)
            for igp in range(0,ngrp):
                F[img,im+igp] = xschi[mat,ig]*xsnf[mat,igp]*volume[ic]

        # end ig in ngrp
    # end j in ncells

    sA0 = sparse.csc_matrix(A0)
    sA2 = sparse.csc_matrix(A2)
    sF  = sparse.csc_matrix(F)

    #plt.spy(A0)
    #plt.show()
    #plt.spy(A2)
    #plt.show()

    if(direct_inv):
        # perform direct matrix inversion
        ## time to do matrix inversion
        print "performing matrix inversion..."
        start_time = time.time()
        A0inv = np.linalg.inv(A0)
        if(method_order == SP3):
            A2inv = np.linalg.inv(A2)
        end_time = time.time()
        print "elapsed time = %.4f" % (end_time-start_time)
    elif(lu_factor):
        # perform LU factorization
        ### time to do LU factorizatoin
        print "performing LU factorization..."
        start_time = time.time()
        #[P0,L0,U0] = scipy.linalg.splu(A0,permute_l=False,overwrite_a=False,check_finite=True)
        #[P2,L2,U2] = scipy.linalg.splu(A2,permute_l=False,overwrite_a=False,check_finite=True)
        luA0 = scipy.sparse.linalg.splu(A0)
        luA2 = scipy.sparse.linalg.splu(A2)
        end_time = time.time()
        print "elapsed time = %.4f" % (end_time-start_time)
        print "done"
    #elif(krylov):
        # do nothing here

    #### solve inverse power iteration CMFD problem

    ####################################
    ######### BEGIN ITERATION ##########
    ####################################

    ###### Outer Power Iteration #######

    t=0
    while(fissdiff > tol and t < 200):
        oldflx[:,:]=scalar_flx[:,:]

        ###################################
        #### nonlineear CMFD iteration ####
        ###################################

        ### no higher order scattering for now
            #for ig in range(0,ngrp):
            #    ### get transport cross sections
            #    xs1p = get_xs1gg(scatmom,matp,ig)
            #    xstrp=xst[matn,ig]-xs1p

            #    xs1n = get_xs1gg(scatmom,matn,ig)
            #    xstrn=xst[matn,ig]-xs1n

            #    dtil[ig,j] = (2.0/3.0)/(xstrn*hn+xstrp*hp)
            #end ig in ngrp
        # end j in ncells+1
       
        #############################
        ##### Setup CMFD Matrix #####
        #############################

        if(krylov):
            rows, cols = M.nonzero()
            for (i,j) in zip(*M.nonzero()):
                M[i,j] = 0.0
        else:
            M.fill(0.0)

        for ic in range(0,ncells):
            im = ic*ngrp
            mat = material_cm[ic]

            ### Set up transport/scattering matrix
            for ig in range(0,ngrp):
                img=im+ig
                aii = 0.0

                xsR0=xst[mat,ig]
                #xsR0=xst[mat,ig]-scatmom[mat][0][ig][ig]
                
                #### loop through all neighbor surfaces, add diffusion term
                for iface in range(0,nfaces):
                    icp = neighbor[iface,ic]

                    dtiln = dtil[0,iface,ig,ic]
                    dhatn = dhat[iface,ig,ic]
                    if(iface == SOUTH or iface == WEST or iface == BOTTOM):
                        dhatn = -dhatn

                    if(icp >= 0):
                        matp = material_cm[icp]
                        hp = h[neighface[iface],icp]
                        ipg = icp*ngrp + ig

                        M[img,ipg] = -dtiln + dhatn
                        aii = aii + dtiln + dhatn
                    else:
                        ### boundary case, do something special
                        aii = aii + dtiln + dhatn

                #### add removal term
                aii = aii + xsR0*volume[ic] 
                M[img,img]=aii

                ### set up off-diagonal scattering 
                gmin=scatbounds[mat][ig][0][0]
                gmax=scatbounds[mat][ig][0][1]
                for igp in range(gmin,gmax):
                    M[img,im+igp]=M[img,im+igp]-scatmom[mat][0][ig][igp]*volume[ic]

                ### fission matrix for CMFD already setup (doesn't change during iteration)
                #for igp in range(0,ngrp):
                #    F[img,im+igp] = xschi[mat,ig]*xsnf[mat,igp]*volume[ic]

            # end ig in ngrp
        # end ic in ncells

        #plt.spy(M)
        #plt.show()
        #plt.spy(F)
        #plt.show()

        #fig = plt.figure()
        #ax = fig.add_subplot(1,1,1)
        #ax.set_aspect('equal')
        #plt.imshow(M, interpolation='none', cmap=plt.cm.ocean)
        #plt.colorbar()
        #plt.show()

        ### print M
        #for ix in range(0,ncells*ngrp):
        #    for iy in range(0,ncells*ngrp):
        #        print " %8.4f" % M[ix,iy],
        #    print " \n",

        #### solve inverse power iteration CMFD problem
        if(direct_inv):
            # direct inversion
            Minv = np.linalg.inv(M)
            invMF = Minv.dot(F)
        elif(lu_factor):
            # LU factorization
            #[PM,LM,UM] = scipy.linalg.lu(M,permute_l=False,overwrite_a=False,check_finite=True)
            luM = scipy.sparse.linalg.splu(M)
        elif(krylov):
            # do nothing here
            sM = sparse.csc_matrix(M)
            ### setup preconditioner
            #MLU = scipy.sparse.linalg.spilu(M)
            #M_x = lambda x: MLU.solve(x)
            #P = scipy.sparse.linalg.LinearOperator((matrix_size,matrix_size), M_x)

        ### calculate keff after transport sweep
        old_fiss_norm = fiss_norm

        keff_cmfd=keff
        lambda_real = 1.0/keff_cmfd
        kdiff=1.0
        wielandt_shift = 0.000
        lamtil = lambda_real - wielandt_shift
        psi_resid=1.0

        tol_cmfd=1E-5
        tol_keff=1E-5
        i=0

        psi = calc_fiss_src_accel(cmflx,material_cm,xsnf,ngrp,ncells,keff_cmfd)

        normphi_init = calc_norm_two(cmflx,ncells*ngrp)
        cmfdstr = "Power Iteration: %i, k_eff = %.5f, k_resid = %.5e, psi_resid = %.5e, normphi = %.4e" % (0,keff_cmfd,tol_keff,tol_cmfd,normphi_init)
        print cmfdstr
        while(psi_resid > tol_cmfd or kdiff > tol_keff):

            #wielandt_shift = 0.5*lambda_real
            #wielandt_shift = 0.0001
            wielandt_shift = 0.0
            lamtil = lambda_real - wielandt_shift

            #   ### apply Wielandt shift
            #Mtil = M
            #for ic in range(0,ncells):
            #    im = ic*ngrp
            #    for ig in range(0,ngrp):
            #        img=im+ig
            #        for igp in range(0,ngrp):
            #            #M[img,im+igp] = M[img,im+igp] - F[img,im+igp]/k_shift
            #            #Mtil[img,im+igp] = M[img,im+igp] - (1.0/k_shift)*F[img,im+igp]
            #            Mtil[img,im+igp] = M[img,im+igp] - wielandt_shift*F[img,im+igp]

            #if(direct_inv):
            #    # direct inversion
            #    Minv = np.linalg.inv(Mtil)
            #    invMF = Minv.dot(F)
            #elif(lu_factor):
            #    # LU factorization
            #    [PM,LM,UM] = scipy.linalg.lu(Mtil,permute_l=False,overwrite_a=False,check_finite=True)
            #elif(krylov):
            #    sM = sparse.csc_matrix(Mtil)
            #    ##sM = M

            if(direct_inv):
                # apply matrix inverse
                #y = invMF.dot(cmflx)/keff_cmfd
                #y = (1.0/keff_cmfd - 1.0/k_shift)*invMF.dot(cmflx)
                y = lamtil*invMF.dot(cmflx)
            elif(lu_factor):
                # solve LU
                #q = sF.dot(cmflx)/keff_cmfd
                #q = (1.0/keff_cmfd - 1.0/k_shift)*sF.dot(cmflx)
                q = lamtil*sF.dot(cmflx)
                #e = np.transpose(PM)
                #u = e.dot(q) 
                #v = scipy.linalg.solve_triangular(LM, u, trans=0, lower=True, unit_diagonal=False, overwrite_b=False, debug=None, check_finite=True)
                #y = scipy.linalg.solve_triangular(UM, v, trans=0, lower=False, unit_diagonal=False, overwrite_b=False, debug=None, check_finite=True)
                y = luM.solve(q)
            elif(krylov):
                # perform Krylov solve
                #b = sF.dot(cmflx)/keff_cmfd
                b = lamtil*sF.dot(cmflx)
                #b = (1.0/keff_cmfd - 1.0/k_shift)*sF.dot(cmflx)*lamtil
                rk = sM.dot(cmflx)-b
                inner_resid = calc_norm_two(rk,ncells*ngrp)
                tmptol = inner_resid*1E-04
                (y, errcode) = scipy.sparse.linalg.gmres(sM, b, x0=cmflx, tol=tmptol,restart=None, maxiter=gmres_max_inner, M=None, callback=None)
                ## perform Krylov solve

            oldk = keff_cmfd

            #### Han Joo's method
            #r1 = y.dot(cmflx)
            #r2 = y.dot(y)
            ##r = r1/r2
            ###keff_cmfd = oldk/r
            ##### shifted
            #gamma = r1/r2
            #lambda_old = lambda_real
            #lambda_real = gamma*lambda_old + (1.0-gamma)*wielandt_shift
            ##keff_cmfd = 1.0/(gamma/oldk + (1.0-gamma)/k_shift)
            ##print "lam_old = %8.3f, lam_shift = %8.3f, lam_new = %8.3f, gamma = %8.3f " % (lambda_old,wielandt_shift,lambda_real,gamma)

            ### MPACT Method
            psi = calc_fiss_src_accel(y,material_cm,xsnf,ngrp,ncells,1.0)
            psid = calc_fiss_src_accel(cmflx,material_cm,xsnf,ngrp,ncells,1.0/lamtil)
            shifted_old = psid+wielandt_shift*psi
            r1 = psi.dot(shifted_old)
            r2 = psi.dot(psi)
            lambda_old = lambda_real
            lambda_real = r1/r2
            lamtil = lambda_real - wielandt_shift

            keff_cmfd = 1.0/lambda_real

            kdiff=abs(keff_cmfd-oldk)

            tmpM = M.dot(y)
            tmpF = F.dot(y)/keff_cmfd
            #tmpF[:] = tmpF[:]*lamtil
            tmpnum = calc_norm_two(tmpM[:] - tmpF[:],ncells*ngrp)
            tmpden = calc_norm_two(tmpF[:],ncells*ngrp)
            phi_resid = tmpnum/tmpden
            #phi_resid=calc_norm_two(tmpresid,ncells*ngrp)
            src_resid=calc_norm_two(tmpF,ncells*ngrp)
            psi_resid=phi_resid/src_resid
            i=i+1
            normphi = calc_norm_two(y,ncells*ngrp)
            cmflx[:] = y[:]

            cmfdstr = "Power Iteration: %i, k_eff = %.5f, k_resid = %.5e, psi_resid = %.5e, normphi = %.4e" % (i,keff_cmfd,kdiff,psi_resid,normphi)
            print cmfdstr

        #### readjust fluxes with CMFD update ####
        ### first, normalize the fluxes 
        for ic in range(0,ncells):
            im=ic*ngrp
            for ig in range(0,ngrp):
                img=im+ig
                scalar_flx[ig,ic]=cmflx[img]
        ### update keff
        keff = keff_cmfd
        #fiss_norm = fiss_norm_new
        fiss_norm = calc_fiss_norm(scalar_flx[:,:],xsnf,volume,material_cm,ncells,ngrp)
        print "fiss_norm = %8.4f " % fiss_norm

        #################################################
        ##### Perform fixed source sweep to update  #####
        ##### estimate of the surface net currents  #####
        #################################################

        ##### Calculate source for fixed-source sweeps #####
        ### sweep over all groups ###
        # calculate the fission source in each cell
        for ic in range(0,ncells):
            mat = material_cm[ic]
            # calculate fission source
            if(xschi[mat,0] > 0.0):
                psi[ic] = calc_psi(oldflx[:,ic],xsnf[mat,:],ig,ngrp)
            else:
            #non fissionable
                psi[ic] = 0.0

        psi[:]=psi[:]/keff

        for isweep in range(0,ninners):

            ### update scattering source for zeroth order ###
            src.fill(0.0)
            for ig in range(0,ngrp):

                # set source in each cell
                for ic in range(0,ncells):
                    mat = material_cm[ic]
                    # calculate scattering source
                    tmpscat = calc_inscat_src(scalar_flx[:,ic],scatmom[mat][0][ig][:],
                        ig,scatbounds[mat,ig,0,0],scatbounds[mat,ig,0,1])
                    ### zeroth order source
                    src[0,ig,ic]=psi[ic]*xschi[mat,ig]+tmpscat

                # end ic in (0,ncells)
            # end ig in (0,ngrp)

            ##### Set fixed source vector #####
            Q0.fill(0.0)
            for ic in range(0,ncells):
                im = ic*ngrp
                mat = material_cm[ic]
                for ig in range(0,ngrp):
                    img = im + ig
                    xsR0 = xst[mat,ig]-scatmom[mat][0][ig][ig]
                    #xsR0 = xst[mat,ig]
                    ### add contribution from second-order
                    Q0[img] = (src[0,ig,ic]+2.0*xsR0*phi[1,ig,ic])*volume[ic]

            if(direct_inv):
                ### use directly inverted matrix
                x = A0inv.dot(Q0)
            elif(lu_factor):
                ### solve LU decomposition
                x = luA0.solve(Q0)
            elif(krylov):
                # perform Krylov solve
                #rk = sA0.dot(x)-Q0
                [x, errcode] = scipy.sparse.linalg.gmres(sA0, Q0, x0=x, tol=1e-06, restart=None, maxiter=gmres_max_inner, M=None, callback=None)
            
            ### place solution back into phi[0,ig,ic]
            for ic in range(0,ncells):
                im = ic*ngrp
                for ig in range(0,ngrp):
                    img = im + ig
                    phi[0,ig,ic] = x[img]
                    #print "ig = %i, ic = %i, phi0 = %8.4f " % (ig,ic,phi[0,ig,ic])

            if(method_order == SP3):
                ##### Calculate source for second-order problem #####
                Q2.fill(0.0)
                for ic in range(0,ncells):
                    im = ic*ngrp
                    mat = material_cm[ic]
                    for ig in range(0,ngrp):
                        img = im + ig
                        xsR0 = xst[mat,ig]-scatmom[mat][0][ig][ig]
                        Q2[img] = 0.4*(xsR0*phi[0,ig,ic] - src[0,ig,ic])*volume[ic]

                ### solve the system
                if(direct_inv):
                    ### use directly inverted matrix
                    x2 = A2inv.dot(Q2)
                elif(lu_factor):
                    ### solve LU decomposition
                    x2 = luA2.solve(Q2)
                elif(krylov):
                    # perform Krylov solve
                    #rk = sA2.dot(x2)-Q2
                    #inner_resid = calc_norm_two(rk,ncells*ngrp)
                    [x2, errcode] = scipy.sparse.linalg.gmres(sA2, Q2, x0=x2, tol=1e-06, restart=None, maxiter=gmres_max_inner, M=None, callback=None)

                ### place solution back into phi[1,ig,ic]
                #print "isweep = %i" % isweep
                for ic in range(0,ncells):
                    im = ic*ngrp
                    for ig in range(0,ngrp):
                        img = im + ig
                        phi[1,ig,ic] = x2[img]
                        scalar_flx[ig,ic] = phi[0,ig,ic] - 2.0*phi[1,ig,ic]

        # end isweep in (0,ninners)


        if(method_order == SP3):
            ### Calculate currents for the following CMFD iteration ###
            for ic in range(0,ncells):
                for iface in range(0,nfaces):
                    icp = neighbor[iface,ic]

                    if(icp >= 0):
                        if(iface == EAST or iface == NORTH or iface == TOP):
                            for ig in range(0,ngrp):
                                surf_cur[iface,ig,ic] = -dtil[0,iface,ig,ic]*(phi[0,ig,icp]-phi[0,ig,ic])
                        elif(iface == WEST or iface == SOUTH or iface == BOTTOM):
                            for ig in range(0,ngrp):
                                surf_cur[iface,ig,ic] = -dtil[0,iface,ig,ic]*(phi[0,ig,ic]-phi[0,ig,icp])
                    else:
                        ### boundary case, do something special
                        if(iface == EAST or iface == NORTH or iface == TOP):
                            for ig in range(0,ngrp):
                                surf_cur[iface,ig,ic] = dtil[0,iface,ig,ic]*phi[0,ig,ic]
                        elif(iface == WEST or iface == SOUTH or iface == BOTTOM):
                            for ig in range(0,ngrp):
                                surf_cur[iface,ig,ic] = -dtil[0,iface,ig,ic]*phi[0,ig,ic]

            # calculate D-hats  
            dhat[:,:,:] = 0.0
            for ic in range(0,ncells):
                for iface in range(0,nfaces):
                    icp = neighbor[iface,ic]
                    if(icp >= 0):
                        for ig in range(0,ngrp):
                            cur=surf_cur[iface,ig,ic]
                            fluxp=scalar_flx[ig,icp]
                            fluxc=scalar_flx[ig,ic]
                            if(iface == EAST or iface == NORTH or iface == TOP):
                                dhat[iface,ig,ic]=(cur+dtil[0,iface,ig,ic]*(fluxp-fluxc))/(fluxp+fluxc)
                            elif(iface == WEST or iface == SOUTH or iface == BOTTOM):
                                dhat[iface,ig,ic]=(cur+dtil[0,iface,ig,ic]*(fluxc-fluxp))/(fluxp+fluxc)
                        # end ig in ngrp
                    else:
                        if(bc[iface] == "reflective"):
                            dhat[iface,:,ic]=0.0
                        elif(bc[iface] == "vacuum"):
                            if(iface == EAST or iface == NORTH or iface == TOP):
                                for ig in range(0,ngrp):
                                    cur=surf_cur[iface,ig,ic]
                                    fluxc=scalar_flx[ig,ic]
                                    dhat[iface,ig,ic]=(cur-dtil[0,iface,ig,ic]*fluxc)/fluxc
                            elif(iface == WEST or iface == SOUTH or iface == BOTTOM):
                                for ig in range(0,ngrp):
                                    cur=surf_cur[iface,ig,ic]
                                    fluxc=scalar_flx[ig,ic]
                                    dhat[iface,ig,ic]=(cur+dtil[0,iface,ig,ic]*fluxc)/fluxc
                                # end ig in ngrp
                # end iface in nfaces
            # end ic in ncells
        elif(method_order == SP1):
            dhat[:,:,:] = 0.0 ### to turn off SP3, effectively

        fissdiff = abs(fiss_norm-old_fiss_norm)
        logstr = "iter: %i, keff = %.5f, kdiff = %.4e, fissdiff = %.4e \n" % (t+1,keff,kdiff,fissdiff)
        logfile.write(logstr)
        print "iter: %i, keff = %.5f, kdiff = %.4e, fissnorm = %.4e, fissdiff = %.4e " % (t+1,keff,kdiff,fiss_norm,fissdiff)
        t=t+1
    # end while fissdiff > tol

    logfile.close()

    ##################
    ### power mesh ###
    ##################

    ### power is printed on 1 cm mesh size. If solution mesh is smaller, 
    ### average powers over the coarse mesh.

    if(nx >= 30):
        nxp = 15
        nyp = 1
        nzp = 15
        dxp = 1.0
        dyp = 1.0
        dzp = 1.0

    else: # too small, don't refine
        nxp = nx/2
        nyp = ny
        nzp = nz/2
        dxp = delta_x
        dyp = delta_y
        dzp = delta_z
    
    nperx = dxp/delta_x   
    npery = dyp/delta_y
    nperz = dzp/delta_z   

    print nperx
    print npery
    print nperz

    nperx = int(nperx)
    npery = int(npery)
    nperz = int(nperz)

    pin_power = np.zeros([nxp,nyp,nzp],dtype=float)
    group_flux = np.zeros([nxp,nyp,nzp,ngrp],dtype=float)
    gflux = np.zeros([ngrp],dtype=float)
    for px in range(0,nxp):
        for py in range(0,nyp):
            for pz in range(0,nzp):
                pinpow = 0.0
                gflux[:] = 0.0
                volsum = 0.0
                ixstt=px*nperx
                iystt=py*npery
                izstt=pz*nperz
                for ix in range(ixstt,ixstt+nperx):
                    for iy in range(iystt,iystt+npery):
                        for iz in range(izstt,izstt+nperz):
                            mat = material_cart[ix,iy,iz]
                            ic = cart_to_cm[ix,iy,iz]
                            volsum=volsum+volume[ic]
                            for ig in range(0,ngrp):
                                pinpow = pinpow + xsf[mat,ig]*scalar_flx[ig,ic]*volume[ic]
                                gflux[ig]+=scalar_flx[ig,ic]*volume[ic]

                pin_power[px,py,pz] = pinpow/volsum
                for ig in range(0,ngrp):
                    group_flux[px,py,pz,ig] = gflux[ig]/volsum

    #########################
    ## normalize pin power ##
    #########################
    powsum = 0.0
    for px in range(0,nxp):
        for py in range(0,nyp):
            for pz in range(0,nzp):
                powsum += pin_power[px,py,pz] 
    normpow = (nxp*nyp*nzp)/powsum
    pin_power[:,:,:] = pin_power[:,:,:]*normpow

    groupsum = 0.0
    for px in range(0,nxp):
        for py in range(0,nyp):
            for pz in range(0,nzp):
                groupsum += group_flux[px,py,pz,0]
    normflux = (nxp*nyp*nzp)/groupsum
    group_flux[:,:,:,:] = group_flux[:,:,:,:]*normflux

    #for ix in range(0,nfuelx):
    #    for iy in range(0,nfuely):
    #        for iz in range(0,nfuelz):
    #            pinpow = 0.0
    #            mat = material_cart[ix,iy,iz]
    #            ic = cart_to_cm[ix,iy,iz]
    #            for ig in range(0,ngrp):
    #                pinpow = pinpow + xsf[mat,ig]*scalar_flx[ig,ic]
    #                group_flux[ix,iy,iz,ig] = scalar_flx[ig,ic]

    #            pin_power[ix,iy,iz] = pinpow

    ###################################
    ### Write Solution to HDF5 File ###
    ###################################

    fuelsize = "%icm" % (int(2.0*delta_x*nfuelx))
    #mesh_delta = "%d" % (delta_x*10.0)
    #meshsize = "%smm" % (mesh_delta)
    h5name = "takeda_%s_%s_h_%s.h5" % (fuelsize,method_name,meshsize)
    hf = h5py.File(h5name,'w')

    g1 = hf.create_group('STATE_0001')
    g1.create_dataset('pin_powers',data=pin_power[:,:,:])
    g1.create_dataset('keff',data=keff)

    g2 = hf.create_group('flux')
    for ig in range(0,ngrp):
        dataname = "group%i" % (ig+1)
        g2.create_dataset(dataname,data=group_flux[:,:,:,ig])

    hf.close()
