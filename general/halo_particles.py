'''
Module containing tools to find and get all particles belonging to
a halo, and to calculate stuff from the full particle information.  
'''

import time
import sys
from numpy import *
import scipy.linalg as linalg
from multiprocessing import Process, Queue, cpu_count
from visnap.general.rockstar import read_rockstar_header
#import pdb

def get_rockstar_halo_particles(rockstar_halo_id, particles_file,
    sanpshot_group, ncpus='all', use_txt_file=0):
    '''
    Get particles belonging to a halo in a Rockstar catalog

    Input:
     rockstar_halo_id - The rockstar catalog halo ID of the halo for which
                        particles are to be found  
 
     particles_file - The rockstar particles files containing the particle
                      information  for the halo with ID = rockstar_halo_id 

     snapshot_group - The snapshot hdf5 group of the IRATE file containing the
                      rockstar halo catalog 
   
     use_txt_file - If set to 1 the particle data will be stored/read in a txt
                    file instead of in a new group within the IRATE file

     ncpus - If 'all' all the available cores in the node will be used, else
             set to the number of cores that you want to use
               

    Output:
     part_data - an numpy array containig the particle data
    '''
    
    if ncpus == 'all':
        ncores = cpu_count()
    else:
        ncores = ncpus
    rockstar_halo_id = int(rockstar_halo_id)
    npart_total = 0 
    
    def find_particles(Nskip,Nread,que,mode=1):
        '''
        Search for particles that belong to the halo

        Kernel function to be called by multiple processes 

        If mode = 0 it only searches for the lines in the txt file that
        contain the particles, with mode = 1  retrieves the data
        '''
        npart = 0
        f = open(particles_file)
        for i,line in enumerate(f):
            if i < Nskip: continue
            if i >= Nskip+Nread: 
                f.close()
                break
            sline = eval(repr(line))
            if '#' not in sline:
                if (int(sline.split()[-1]) == rockstar_halo_id):
                    if mode:
                        x,y,z,vx,vy,vz,particle_id,ai_id,i_id,e_id = sline.split()
                        if npart == 0:
                            part_data = array([x,y,z,vx,vy,vz,particle_id,ai_id,i_id,e_id]).reshape(1,10)
                        else:
                            this_part_data = array([x,y,z,vx,vy,vz,particle_id,ai_id,i_id,e_id]).reshape(1,10)
                            part_data = concatenate((part_data,this_part_data))
                        npart += 1
                    else:
                        if npart == 0: part_data = (i,-1)
                        else: part_data = (part_data[0],i)
                        npart += 1
               
        if npart > 0:
            que.put(part_data)
            return      
        else: 
            que.put(-1)
            return    

    if use_txt_file:
        outfile = particles_file.replace('particles','halo'+str(rockstar_halo_id)+'_particles')    
        if os.path.exists(outfile):
            print "Found file %s containing particle data for halo %d" % (outfile, rockstar_halo_id)
            part_data = loadtxt(outfile)
            npart_total = part_data[:,0].size
            print "Found %d particles belonging to halo %d" % (npart_total,rockstar_halo_id)
            return part_data
        else:
            fin = open(particles_file,'r')
            fout = open(outfile,'w')
            # write header to outfile
            for i,line in enumerate(fin):
                sline = eval(repr(line))
                if '#Halo table begins here' not in sline: fout.write(sline)
                else: 
                    fout.write('#Particle table begins here:\n')
                    break
            fin.close()
    else: 
        #using IRATE file to read/write particle data
        snap = sanpshot_group
        if 'halo'+str(rockstar_halo_id)+'_particles' in snap.keys():
            print "Found data set containing particle data for halo %d in group %s" % (rockstar_halo_id, snap.name)
            pd_group = snap['halo'+str(rockstar_halo_id)+'_particles']
            part_data = array([pd_group['Position'][...][:,0],
                               pd_group['Position'][...][:,1],
                               pd_group['Position'][...][:,2],
                               pd_group['Velocity'][...][:,0],
                               pd_group['Velocity'][...][:,1],
                               pd_group['Velocity'][...][:,2],
                               pd_group['particle_id'][...],
                               pd_group['assigned_internal_haloid'][...],
                               pd_group['internal_haloid'][...],
                               pd_group['external_haloid'][...]])
            part_data = part_data.transpose()
            npart_total = part_data[:,0].size
            print "Found %d particles belonging to halo %d" % (npart_total,rockstar_halo_id)
            return part_data
        else:
            #write comments and cosmology as attributes, also create group and prepare labels/units stuff 
            comments, cosmo, unitnames, unitcgs  = read_rockstar_header(particles_file)
            pd_group = snap.create_group('halo'+str(rockstar_halo_id)+'_particles')
            labels = ['Position','Velocity','particle_id',
                      'assigned_internal_haloid','internal_haloid',
                      'external_haloid']
            units_names = [unitnames['Position'],unitnames['Velocity'],None,None,None,None]
            units_cgs = [unitcgs['Position'],unitcgs['Velocity'],None,None,None,None]
            attrs = comments+cosmo
            for a in attrs:
                try:
                    a[1] = float(a[1])
                except ValueError:
                    pass
            for i in range(len(attrs)):
                pd_group.attrs[attrs[i][0]] = attrs[i][1]


        print "Looking for particles in file %s that belong to halo %d" % (particles_file,rockstar_halo_id)
        time.sleep(3)
        fin = open(particles_file,'r')
        nlines = sum(1 for line in fin)
        nlines_pc = nlines/ncores
        nleft = nlines%ncores
        fin.close()

        que = Queue()
        jobs = []
        for i in range(ncores):
            Nskip = nlines_pc*i
            if i < ncores-1: Nread = nlines_pc
            else: Nread = nlines_pc+nleft
            #first we run in mode = 0 to find the lines in the file that have the particles we want 
            thisJob = Process(target = find_particles, args = (Nskip,Nread,que,0))
            jobs.append(thisJob)

        for job in jobs: job.start()
        while que.qsize() < len(jobs): time.sleep(5)
        results = [que.get() for i in range(que.qsize())]
        for job in jobs: job.join()

        first_line = Inf
        last_line = 0
        for result in results:
            if type(result) != type(-1):
                first_line = min(first_line,result[0])
                last_line = max(last_line,result[1])
                
        if last_line == 0:
            print 'No particles found in %s belonging to halo %d' % (particles_file,rockstar_halo_id)
            if use_txt_file:
                fout.close()
                os.system('rm '+ outfile)
            sys.exit()
        else:
            print 'Particles found! reading %d lines of data with %d cores'% (last_line-first_line,ncores)
            nlines = last_line-first_line
            nlines_pc = nlines/ncores
            nleft = nlines%ncores

            que = Queue()
            jobs = []
            for i in range(ncores):
                Nskip = first_line + nlines_pc*i
                if i < ncores-1: Nread = nlines_pc
                else: Nread = nlines_pc+nleft
                #Now we get the particle data reading only the lines where the particles we want are
                thisJob = Process(target = find_particles, args = (Nskip,Nread,que)) 
                jobs.append(thisJob)
        
            for job in jobs: job.start()
            while que.qsize() < len(jobs): time.sleep(5)
            results = [que.get() for i in range(que.qsize())]
            for job in jobs: job.join()

            FirstResultFlag = 1
            for result in results:
                if type(result) != type(-1):
                    if  FirstResultFlag:
                        part_data = result
                        FirstResultFlag = 0
                    else: part_data = concatenate((part_data,result))
            
            npart_total = part_data[:,0].size
            print "Found %d particles belonging to halo %d" % (npart_total,rockstar_halo_id)
            
            if use_txt_file:
                print "Writing particle data to file %s"%(outfile)        
                for i in range(npart_total):
                    fout.write('%s %s %s %s %s %s %s %s %s %s\n' % (part_data[i,0],part_data[i,1],part_data[i,2],part_data[i,3],part_data[i,4],part_data[i,5],part_data[i,6],part_data[i,7],part_data[i,8],part_data[i,9]))
                fout.close()
            else:
                print "Saving particle data in the open IRATE file being analyzed"
                print "Writing particle data under %s" % (pd_group.name)
                pos, vel = array([part_data[:,0],part_data[:,1],part_data[:,2]]), array([part_data[:,3],part_data[:,4],part_data[:,5]])
                pos, vel = pos.transpose(), vel.transpose()
                data_sets = [pos,vel,part_data[:,6],part_data[:,7],part_data[:,8],part_data[:,9]]
                for i in range(len(labels)):
                    pd_group.create_dataset(labels[i],data=data_sets[i])
                    if units_names[i] != None:
                        pd_group[labels[i]].attrs['unitname'] = units_names[i]
                    if units_cgs[i] != None:
                        pd_group[labels[i]].attrs['unitcgs'] = units_cgs[i]

            return part_data



def calculate_profiles(part_data, hcenter, hvelocity, nbins=15):
    '''
    Generate radial profiles from full particle data

    Input:

     part_data - particle data array [x,y,z,vx,vy,vz]

     hcenter - the halo center

     hvelocity - the halo bulk velocity

     nbins - number of bins in r

    Output:

     rmid - the mid point r for that r bin

     Ninshell - number of particles in bin

     Nenclosed - number of particles with r<rmid

     dens - number density for that bin 

     avgDens - average number density inside rmid (eclosed density)

     vdisp - total velocity dispersion profile
    '''
    # ignore divide by zero and invalid operation warnings
    seterr(divide='ignore', invalid='ignore') 

    part_data,hcenter,hvelocity = part_data.astype(float),hcenter.astype(float),hvelocity.astype(float)
    x,y,z = part_data[:,0]-hcenter[0],part_data[:,1]-hcenter[1],part_data[:,2]-hcenter[2]  
    vx,vy,vz = part_data[:,3]-hvelocity[0],part_data[:,4]-hvelocity[1],part_data[:,5]-hvelocity[2]
    r = sqrt(x*x + y*y + z*z) 
    
    radbins = logspace(log10(r[r != 0].min()),log10(r.max()),num=nbins)
    rmid =  zeros(len(radbins)-1,dtype='float')
    Ninshell = zeros(len(radbins)-1)
    Nenclosed = zeros(len(radbins)-1)
    dens = zeros(len(radbins)-1,dtype='float')     
    avgDens = zeros(len(radbins)-1,dtype='float')     
    vdisp = zeros(len(radbins)-1,dtype='float')    
    
    for i in range(len(radbins)-1):      
        lowr = radbins[i]
        highr = radbins[i+1]
        rmid[i] = (lowr+highr)/2.0
        inshell =  (r > lowr) & (r < highr)
        
        # Differential density 
        Ninshell[i] = r[inshell].size     
        vol = (highr**3 - lowr**3)*4.0*pi/3.0       
        dens[i] = Ninshell[i]/vol  # This is a number density      
       
        # Average (enclosed) density
        Nenclosed[i] = r[r < rmid[i]].size
        avgDens[i] = Nenclosed[i]/(4.0*pi*rmid[i]**3/3.0)
       
        # Total Velocity Dispersion
        meanVx = vx[inshell].mean()
        meanV2x = (vx[inshell]*vx[inshell]).mean()
        vxdisp = sqrt(meanV2x - meanVx*meanVx)

        meanVy = vy[inshell].mean()
        meanV2y = (vy[inshell]*vy[inshell]).mean()
        vydisp = sqrt(meanV2y - meanVy*meanVy)

        meanVz = vz[inshell].mean()
        meanV2z =  (vz[inshell]*vz[inshell]).mean()
        vzdisp = sqrt(meanV2z - meanVz*meanVz)

        vdisp[i] = sqrt(vxdisp**2 + vydisp**2 + vzdisp**2)
        
    return rmid, Ninshell, Nenclosed, dens, avgDens, vdisp



def calculate_2dprofiles(part_data, hcenter, hvelocity, projection_axis=(1,0,0), nbins=15):
    '''
    Generate radial 2d projection profiles from full particle data

    Input:

     part_data - particle data array [x,y,z,vx,vy,vz]

     hcenter - the halo center
     
     hvelocity - the halo bulk velocity

     projection_axis - The direction of the projection axis in a vector form
                       i.e. (x, y, z)
     
     nbins - number of bins in r

    Output:

     Rmid - the mid point R for that R bin

     Ninshell - number of particles in bin

     Nenclosed - number of particles with R<Rmid

     Sdens - number surface density for that bin 

     avgSdens - average number surface density inside Rmid (eclosed number density)

     vlosDisp - Line of sight velocity dispersion profile

     vtanDisp - Tangential velocity dispersion profile
    '''
    # ignore divide by zero and invalid operation warnings
    seterr(divide='ignore', invalid='ignore') 

    part_data,hcenter,hvelocity = part_data.astype(float),hcenter.astype(float),hvelocity.astype(float)
    r_vec = part_data[:,0:3]-hcenter
    v_vec = part_data[:,3:6]-hvelocity
    r_mag = array([sqrt(dot(thisr_vec,thisr_vec)) for thisr_vec in r_vec])
    v_mag = array([sqrt(dot(thisv_vec,thisv_vec)) for thisv_vec in v_vec])

    eproj = projection_axis/sqrt(dot(projection_axis,projection_axis))  # normalize 
    r_dot_eproj = array([dot(thisr_vec,eproj) for thisr_vec in r_vec])
    vlos = array([dot(thisv_vec,eproj) for thisv_vec in v_vec])
    vtan = sqrt(v_mag*v_mag - vlos*vlos)
    R = sqrt(r_mag*r_mag - r_dot_eproj*r_dot_eproj)   
    
    radbins = logspace(log10(R[R != 0].min()),log10(R.max()),num=nbins)
    Rmid =  zeros(len(radbins)-1,dtype='float')
    Ninshell = zeros(len(radbins)-1)
    Nenclosed = zeros(len(radbins)-1)
    Sdens = zeros(len(radbins)-1,dtype='float')     
    avgSdens = zeros(len(radbins)-1,dtype='float')     
    losVdisp = zeros(len(radbins)-1,dtype='float')    
    
    for i in range(len(radbins)-1):      
        lowR = radbins[i]
        highR = radbins[i+1]
        Rmid[i] = (lowR+highR)/2.0
        inshell =  (R > lowR) & (R < highR)
        
        # Differential density 
        Ninshell[i] = R[inshell].size     
        area = pi*(highR**2 - lowR**2)       
        Sdens[i] = Ninshell[i]/area  # This is a number density      
       
        # Average (enclosed) density
        Nenclosed[i] = R[R < Rmid[i]].size
        avgSdens[i] = Nenclosed[i]/(pi*Rmid[i]**2)
       
        # Velocity Dispersion
        meanVlos = vlos[inshell].mean()
        meanVlos2 = (vlos[inshell]*vlos[inshell]).mean()
        vlosDisp = sqrt(meanVlos2 - meanVlos*meanVlos)
        
        meanVtan = vtan[inshell].mean()
        meanVtan2 = (vtan[inshell]*vtan[inshell]).mean()
        vtanDisp = sqrt(meanVtan2 - meanVtan*meanVtan)

    return Rmid, Ninshell, Nenclosed, Sdens, avgSdens, vlosDisp, vtanDisp



def axis_ratios(arr_in, rad, shell=False, axes_out=False, fix_volume=True, quiet=False):
    """Compute axis ratios iteratively.

    WORK IN PROGRESS -- converting from gadget_profile.pro.  
    Needs to be checked.

    May want to add capability to compute axis ratios at multiple
    different radii in one function call (w/ a loop?) 

    See Dubinski & Carlberg (1991) for details

    INPUTS:
    arr_in: array of particle positions, assumed to be centered
    rad: radius to compute axis ratios.  If computing axis ratios in
    shells, rad should be a two-element list / array, with rad[0] =
    inner radius of shell and rad[1] = outer radius of shell.
    Otherwise, rad should be a real number equal to the radius within
    which to compute the axis ratios


    OPTIONAL INPUTS:
    shell=False: by default, compute cumulative axis ratio
    (i.e. for all particles within radius r).  If shell=True, compute
    axis ratio in an ellipsoidal shell instead.
    axes_out=False:  if True, also return principal axes (in ascending order)
    fix_volume=True: keep the volume of the ellipsoid constant while iterating.
    If false, keep the semi-major axis equal to the initial (spherical) search
    radius. This will result in a smaller effective volume.
    quiet=False: if set to true, suppress information that is printed
    on each iteration
    """
    def calc_inertia(arr_in, axrat_in):
	"""calculate the modified moment of inertia tensor and get its
	eigenvalues and eigenvalues"""
	tensor=zeros([3,3])
	# given initial values for primary axes, compute ellipsoidal
	# radius of each particle
	rp2=(arr_in[:,0]**2 +
             arr_in[:,1]**2/axrat_in[0]**2 + 
             arr_in[:,2]**2/axrat_in[1]**2)

	# construct the moment of inerial tensor to be diagonalized:
	for i in range(3):
	    for j in range(3):
		tensor[i,j]=(arr_in[:,i]*arr_in[:,j]/rp2).sum()
	
	evecs=linalg.eig(tensor)
    
	# return a tuple with eigenvalues and eigenvectors; 
	# eigenvalues=evecs[0], eigenvectors=evecs[1]
	return evecs

    cnt=0
    # initial guess for principal axes:
    evs0=array([[1e0,0e0,0e0],[0e0, 1e0, 0e0], [0e0, 0e0, 1e0]])
    # initial guess for axis ratios:
    axes=array([1e0,1e0])
    avdiff=1e0

    rad=asarray(rad)
    while (avdiff > 0.01) & (cnt < 100):
	axtemp=axes.copy()
	# compute ellipsoidal radius of each particle:
	dist2=(arr_in[:,0]**2 +
	       arr_in[:,1]**2/axes[0]**2 +
	       arr_in[:,2]**2/axes[1]**2)**0.5

	if shell == True:
	    # find locations of particles between rad[0] and rad[1]
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**3/axes[0]/axes[1])**(1e0/3e0)
	    locs=((dist2 < r_ell[1]) & 
                  (dist2 > r_ell[0])).nonzero()[0]
	else:
	    # find locations of particles with r < rad
	    locs=(dist2 < rad).nonzero()
	    # current version keeps volume same as original spherical
	    # volume rather than keeping principal axis same as
	    # initial radius
	    # compute ellipsoidal radius = (a*b*c)^(1/3)
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**3/axes[0]/axes[1])**(1e0/3e0)
	# get eigenvectors and eigenvalues.  Note: the inertia tensor
	# is a real symmetric matrix, so the eigenvalues should be real.
	axrat=calc_inertia(arr_in[locs], axtemp)
	if abs(imag(axrat[0])).max() > 0.:
	    #raise ValueError('Error: eigenvalues are not all real!')
            #Annika edit so that this doesn't totally crap out if I have to loop over it.
            print 'Error: eigenvalues are not all real!'
            if axes_out:
                return [-1., -1.], [evs0[:,inds[0]], evs0[:,inds[1]], evs0[:,inds[2]]]
            else:
                return [-1., -1.]
            

	evals=real(axrat[0])
	evecs=axrat[1]
	# get axis ratios from eigenvalues:
	axes=sqrt([evals[1]/evals[0], 
                      evals[2]/evals[0]])
	# rotate particles (and previous eigenvector array) into new
	# basis:
	arr_in=dot(arr_in, evecs)
	evs0=dot(evs0, evecs)
	
	# compute difference between current and previous iteration
	# for axis ratios and diagonal of eigenvector matrix
	avd0=abs((axtemp[0]-axes[0])/(axtemp[0]+axes[0]))
	avd1=abs((axtemp[1]-axes[1])/(axtemp[1]+axes[1]))

	# used to check how close most recent rotation is to zero
	# rotation (i.e. how close evecs is to the identity matrix)
	avd2=3e0-abs(evecs[0,0])-abs(evecs[1,1]) - abs(evecs[2,2])
	if quiet == False:
	    # print axis ratios relative to major axis
	    print 'axis ratios: ', (sort(evals/evals.max()))**0.5
            print 'deviations from previous iteration: ' + \
		'%.*e, %.*e, %.*e' % (4, avd0, 4, avd1, 4, avd2)
	    print 'number of particles in shell / sphere: ', size(locs)
	avdiff=max(avd0, avd1, avd2)
	cnt+=1

    # normalize eigenvalues to largest eigenvalue
    evals=evals/evals.max()
    inds=evals.argsort()
    final_axrats=evals[inds[:2]]**0.5
    print 'number of iterations: ', cnt
    print 'number of particles in shell / sphere: ', size(locs), '\n'
    print 'normalized eigenvalues (sorted in ascending order) ' + \
	'and corresponding unit eigenvectors are:'
    print '%.*f' % (5, evals[inds[0]]), evs0[:,inds[0]]
    print '%.*f' % (5, evals[inds[1]]), evs0[:,inds[1]]
    print '%.*f' % (5, evals[inds[2]]), evs0[:,inds[2]]
    print 'axis ratios: ', final_axrats
    if cnt == 100:
        print 'Failed to converge after 100 iterations'
    if axes_out:
        return final_axrats, [evs0[:,inds[0]], evs0[:,inds[1]], evs0[:,inds[2]]]
    else:
        return final_axrats


"""Annika edit for 2d shapes"""
"""Annika Peter 12/2/2011"""
def axis_ratios2D(arr_in, rad, shell=False, axes_out=False, fix_volume=True, quiet=False):
    """Compute axis ratios iteratively.

    WORK IN PROGRESS -- converting from gadget_profile.pro.  
    Needs to be checked.

    May want to add capability to compute axis ratios at multiple
    different radii in one function call (w/ a loop?) 

    See Dubinski & Carlberg (1991) for details

    INPUTS:
    arr_in: array of particle positions, assumed to be centered
    rad: radius to compute axis ratios.  If computing axis ratios in
    shells, rad should be a two-element list / array, with rad[0] =
    inner radius of shell and rad[1] = outer radius of shell.
    Otherwise, rad should be a real number equal to the radius within
    which to compute the axis ratios


    OPTIONAL INPUTS:
    shell=False: by default, compute cumulative axis ratio
    (i.e. for all particles within radius r).  If shell=True, compute
    axis ratio in an ellipsoidal shell instead.
    axes_out=False:  if True, also return principal axes (in ascending order)
    fix_volume=True: keep the volume of the ellipsoid constant while iterating.
    If false, keep the semi-major axis equal to the initial (spherical) search
    radius. This will result in a smaller effective volume.
    quiet=False: if set to true, suppress information that is printed
    on each iteration
    """
    def calc_inertia(arr_in, axrat_in):
	"""calculate the modified moment of inertia tensor and get its
	eigenvalues and eigenvalues"""
	tensor=zeros([2,2])
	# given initial values for primary axes, compute ellipsoidal
	# radius of each particle
	rp2=(arr_in[:,0]**2 +
             arr_in[:,1]**2/axrat_in[0]**2)

	# construct the moment of inerial tensor to be diagonalized:
	for i in range(2):
	    for j in range(2):
		tensor[i,j]=(arr_in[:,i]*arr_in[:,j]/rp2).sum()
	
	evecs=linalg.eig(tensor)
    
	# return a tuple with eigenvalues and eigenvectors; 
	# eigenvalues=evecs[0], eigenvectors=evecs[1]
	return evecs

    cnt=0
    # initial guess for principal axes:
    evs0=array([[1e0,0e0],[0e0, 1e0]])
    # initial guess for axis ratios:
    axes=array([1e0])
    avdiff=1e0

    rad=asarray(rad)
    while (avdiff > 0.01) & (cnt < 100):
	axtemp=axes.copy()
	# compute ellipsoidal radius of each particle:
	dist2=(arr_in[:,0]**2 +
	       arr_in[:,1]**2/axes[0]**2)**0.5

	if shell == True:
	    # find locations of particles between rad[0] and rad[1]
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**2/axes[0])**(1e0/2e0)
	    locs=((dist2 < r_ell[1]) & 
                  (dist2 > r_ell[0])).nonzero()[0]
	else:
	    # find locations of particles with r < rad
	    locs=(dist2 < rad).nonzero()
	    # current version keeps volume same as original spherical
	    # volume rather than keeping principal axis same as
	    # initial radius
	    # compute ellipsoidal radius = (a*b*c)^(1/3)
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**2/axes[0])**(1e0/2e0)
	# get eigenvectors and eigenvalues.  Note: the inertia tensor
	# is a real symmetric matrix, so the eigenvalues should be real.
	axrat=calc_inertia(arr_in[locs], axtemp)
	if abs(imag(axrat[0])).max() > 0.:
	    #raise ValueError('Error: eigenvalues are not all real!')
            #Annika edit so that this doesn't totally crap out if I have to loop over it.
            print 'Error: eigenvalues are not all real!'
            if axes_out:
                return [-1.], [evs0[:,inds[0]], evs0[:,inds[1]]]
            else:
                return [-1.]
            

	evals=real(axrat[0])
	evecs=axrat[1]
	# get axis ratios from eigenvalues:
	axes=sqrt([evals[1]/evals[0]])
	# rotate particles (and previous eigenvector array) into new
	# basis:
	arr_in=dot(arr_in, evecs)
	evs0=dot(evs0, evecs)
	
	# compute difference between current and previous iteration
	# for axis ratios and diagonal of eigenvector matrix
	avd0=abs((axtemp[0]-axes[0])/(axtemp[0]+axes[0]))

	# used to check how close most recent rotation is to zero
	# rotation (i.e. how close evecs is to the identity matrix)
	avd2=2e0-abs(evecs[0,0])-abs(evecs[1,1])
	if quiet == False:
	    # print axis ratios relative to major axis
	    print 'axis ratios: ', (sort(evals/evals.max()))**0.5
            print 'deviations from previous iteration: ' + \
		'%.*e, %.*e' % (4, avd0, 4, avd2)
	    print 'number of particles in shell / sphere: ', size(locs)
	avdiff=max(avd0, avd2)
	cnt+=1

    # normalize eigenvalues to largest eigenvalue
    evals=evals/evals.max()
    inds=evals.argsort()
    final_axrats=evals[inds[:2]]**0.5
    print 'number of iterations: ', cnt
    print 'number of particles in shell / sphere: ', size(locs), '\n'
    print 'normalized eigenvalues (sorted in ascending order) ' + \
	'and corresponding unit eigenvectors are:'
    print '%.*f' % (5, evals[inds[0]]), evs0[:,inds[0]]
    print '%.*f' % (5, evals[inds[1]]), evs0[:,inds[1]]
    print 'axis ratios: ', final_axrats
    if cnt == 100:
        print 'Failed to converge after 100 iterations'
    if axes_out:
        return final_axrats, [evs0[:,inds[0]], evs0[:,inds[1]]]
    else:
        return final_axrats
