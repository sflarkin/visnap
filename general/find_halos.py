'''Module containing tools for finding specfic halos in a simulation'''

import sys
import h5py
from numpy import * 
from visnap.general import halo
#import pdb

def find_zoom_halo(halo_list, zoom_id, irate_file, snapshot, catalog,  
                   prop='Mvir', prop_value=None, N_mostparticles=5):
    '''
    Find the halo in the given catalog that most closely matches the given property
    or the properties of the halo with zoom_id given in the halo_list file. This is 
    designed to find Zoom halos in the hi-res region of a simulation, thus a selection
    is done first consisting of the N_mostparticles halos with the most number
    of particles in the catalog, the best match within this selction will be
    returned.  

    Input: 
    
     halo_list - A file containing a list of zoom halos with their zoom IDs and
                 halo properties (if prop_value is not None halo_list will be ignored)
    
     zoom_id - A zoom halo ID (as set in the IRATE filename for a zoom
               simulation). Note that this is normally different to the halo ID
               in the catalog (if prop_value is not None zoom_id will be ignored)
     
     irate_file - The name of the IRATE file containing the catalog

     snapshot - The snapshot from which the catalog was generated.
                Can be a string with the h5py key or an integer
                with the snapshot number          

     catalog - The AHF or Rockstar catalog. Can be a string with 
               the h5py key or simply "AHF"/"Rockstar", if given as
               the later the first "AHF"/"Rockstar" catalog found will
               be used

     prop - property to use for the matching, can be "Mvir" (default), "Vmax",
            "Rvir", "Vel", "Spin" or "all"

     prop_value - If given, halo_list and zoom_id will be ignored, and  
                  prop_value will be used as the value for the halo property
                  to match
     
     N_mostparticles - Number of halos in the initial selection of halos with
                       the most number of particles

    Output:

     halo_object - The halo object of the found halo as created by
                   visnap.halo.new_halo() 
              
     dist - A distance measure of how closely the properties were matched
    '''
   
    # Get data from the zoom list
    if prop_value == None:
        zoom_list = genfromtxt(halo_list)
        ids = zoom_list.T[0]
        arg = argwhere(ids == int(zoom_id))[0,0]
        hid_l, Vx_l, Vy_l, Vz_l, Rvir_l, Mvir_l, Vmax_l, Spin_l = zoom_list[arg]

    # open snapshot
    irate = h5py.File(irate_file,'r')
    if type(snapshot) == int: snapshot = 'Snapshot%.5d' % snapshot
    snap = irate[snapshot]   
    
    # open catalog
    C = catalog
    if (catalog == 'AHF') or (catalog == 'Rockstar'):
        for key in snap.keys():
            if catalog in key:
                C = snap[key]
                catalog = key
                break
            if key == snap.keys()[-1]:
                errmsg = 'No %s catalog found' % catalog
                raise ValueError(errmsg)
    else: 
        C = snap[catalog]

    # Get the properties of the N_mostparticles halos in the catalog    
    npart = C['Num_p'][...]
    args = npart.argsort()[-N_mostparticles:]
    Vmax, Mvir, Rvir = [C['Vmax'][...][args], C['Mvir'][...][args],
                        C['Rvir'][...][args]]
    Vx, Vy, Vz = [C['Velocity'][...][:,0][args], C['Velocity'][...][:,1][args],
                  C['Velocity'][...][:,2][args]]
    try:
        Spin = C['lambda'][...][args]
        spin_flag = 1
    except KeyError:
        try:
            Spin = C['Spin'][...][args]
            spin_flag = 1
        except KeyError:
            print "No spin information found in this catalog"
            spin_flag = 0
            if prop == 'Spin':
                print "Try using a different propertie to match halo"
                sys.exit()
            
    try:
        ID =  C['ID'][...][args]
        id_flag = 1
    except KeyError:
        id_flag = 0
   
    # Compute distances    
      
    if prop == 'Mvir':
        if prop_value == None: MvirDist = (Mvir-Mvir_l)/Mvir_l
        else:  MvirDist = (Mvir-prop_value)/prop_value
        dist = sqrt(MvirDist**2)
    elif prop == 'Vmax':
        if prop_value == None: VmaxDist = (Vmax-Vmax_l)/Vmax_l
        else:  VmaxDist = (Vmax-prop_value)/prop_value
        dist = sqrt(VmaxDist**2)
    elif prop == 'Rvir': 
        if prop_value == None: RvirDist = (Rvir-Rvir_l)/Rvir_l
        else: RvirDist = (Rvir-prop_value)/prop_value
        dist = sqrt(RvirDist**2)
    elif prop == 'Vel':
        if prop_value == None: 
            VxDist, VyDist, VzDist = (Vx-Vx_l)/Vx_l,(Vy-Vy_l)/Vy_l,(Vz-Vz_l)/Vz_l
        else:
            VxDist, VyDist, VzDist = (Vx-prop_value[0])/prop_value[0],(Vy-prop_value[1])/prop_value[1],(Vz-prop_value[2])/prop_value[2]
        VDist = sqrt(VxDist**2 + VyDist**2 + VzDist**2)    
        dist = VDist
    elif prop == 'Spin':
        if prop_value == None: SpinDist = (Spin-Spin_l)/Spin_l
        else:  SpinDist = (Spin-prop_value)/prop_value
        dist = sqrt(SpinDist**2)
    elif prop == 'all':  
        if prop_value != None:
            print "A prop_value was given, thus you can not use prop='all'"
            sys.exit()
        if spin_flag:
            dist = sqrt( VmaxDist**2 + MvirDist**2 + RvirDist**2 +
                            SpinDist**2 + VxDist**2 + VyDist**2 + VzDist**2 )
        else:
            dist = sqrt( VmaxDist**2 + MvirDist**2 + RvirDist**2 +
                            VxDist**2 + VyDist**2 + VzDist**2 )
    else:
        print 'Unrecognized property %s to match'%prop
        sys.exit()
    
    # Get halo_arg and halo_id    
    arg =  argwhere(dist == dist.min())[0,0]
    halo_arg = argwhere(C['Vmax'][...] == Vmax[arg])[0,0]
    if id_flag:
        hal_id = C['Vmax'][...][halo_arg]
    
    # Create halo_object    
    halo_object = halo.new_halo(irate_file, snapshot, catalog,
                                halo_arg=halo_arg)
    #Done
    irate.close()
    return halo_object, dist.min()


def find_mostpart_halo(irate_file, snapshot, catalog, nhalos=1):
    '''
    Find the halos with the largest number of particles
 
    Input:
  
     irate_file - The name of the IRATE file containing the catalog

     snapshot - The snapshot from which the catalog was generated.
                Can be a string with the h5py key or an integer
                with the snapshot number. If snapshot=None then the
                file will be treated as an old IRATE.

     catalog - The AHF or Rockstar catalog. Can be a string with 
               the h5py key or simply "AHF"/"Rockstar", if given as
               the later the first AHF/Rockstar catalog found will
               be used

     nhalos - If nhalos > 1 will return a list with the nhalos 
              that have the largest number of particles

    Output:

     halo_objects - If nhalos > 1: A list of halo objects as created by
                    visnap.halo.new_halo(). Else: A single halo object. 
    ''' 

    # open snapshot
    irate = h5py.File(irate_file,'r')
    if snapshot:
        if type(snapshot) == int: snapshot = 'Snapshot%.5d' % snapshot
        snap = irate[snapshot]
        # open catalog
        if (catalog == 'AHF') or (catalog == 'Rockstar'):
            for key in snap.keys():
                if catalog in key:
                    C = snap[key]
                    catalog = key
                    break
                if key == snap.keys()[-1]:
                    errmsg = 'No %s catalog found' % catalog
                    raise ValueError(errmsg)
        else: 
            C = snap[catalog]
    else:
        C = irate[catalog]

    # Find halos with the largest number of particles    
    npart = C['Num_p'][...]
    halo_args = npart.argsort()[-nhalos:][::-1]

    halos = []
    # Create halo_objects
    for arg in halo_args:
        thishalo = halo.new_halo(irate_file, snapshot, catalog,
                                 halo_arg=arg)
        halos.append(thishalo)
        
    # Done
    irate.close()
    if len(halos)==1:
        return halos[0]
    else:
        return halos


def find_subhalos(halo_object, vcut_factor=0.05, rcut_factor=1.0, min_npart=100):
    '''
    Find the subhalos of the given halo 

    Input
 
     halo_object - A halo object as returned by visnap.general.halo.new_halo()

     vcut_factor - The lowest subhalo vmax to consider in units of the host Vmax

     rcut_factor - Only subhalos with r < rcut_factor*host_Rvir will be returned            
  
     min_npart - Only subhalos with at least min_npart number of particles
                 will be included 
    
    Output

     subhalo_list - A list of halo objects with all the found subhalos
    '''

    # open snapshot
    host = halo_object
    irate = h5py.File(host.irate_file,'r')
    snap = irate[host.snapshot]

    # open catalog
    C = snap[host.catalog]

    try: 
        hosts = C['Hosts'][...]
    except KeyError:
        try:
            hosts = C['hostHalo'][...]
        except KeyError:
            print "No host/subhalo identifyer found in this IRATE catalog. "\
                "Now generating, unbound subhalos are NOT being removed" 
            from visnap.general import subhalo_identify
            subhalo_identify(C)
            hosts = C['Hosts'][...]

    npart = C['Num_p'][...]  
    vmax = C['Vmax'][...] 
    center  = C['Center'][...]
    center_units = C['Center'].attrs['unitname']
    host_center = host.props['Center']
    host_velocity = host.props['Velocity']
    host_rvir, host_vmax = host.props['Rvir'], host.props['Vmax']
    
    
    R = sqrt((center[:,0] - host_center[0])**2 +
             (center[:,1] - host_center[1])**2 +
             (center[:,2] - host_center[2])**2)
    if 'Mpc' in center_units:
        R = R*1000

    print "Selecting all subhalos with more than %d particles within %gRvir "\
        "of the host and with vmax > %g km/s" % (min_npart, rcut_factor, vcut_factor*host_vmax)
    sub_cut = argwhere((hosts != -1) & (npart > min_npart) & 
                       (vmax > vcut_factor*host_vmax) & 
                       (R < rcut_factor*host_rvir) & (R > 0))[:,0] 
    
    subhalo_list = []
    print 'Generating subhalo list'
    for arg in sub_cut:
        subhalo_object = halo.new_halo(host.irate_file, host.snapshot,
                                       host.catalog, halo_arg=arg)
        subhalo_list.append(subhalo_object)

    irate.close()    
    return subhalo_list
                                       
        
        
    
