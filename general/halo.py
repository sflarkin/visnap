'''Module containing tools for generating halo objects'''

import sys
import h5py
from numpy import *
import visnap
import visnap.general.irate_file_mod as irate_file_mod
#import pdb

class new_halo:
    '''Create a new halo object'''

    def __init__(self, irate_file, snapshot, catalog, halo_id=None, 
                 halo_arg=None):
        '''
        Initialize new halo object with basic attributes and properties

        irate_file - The name of the IRATE file containing the catalog

        snapshot - The snapshot from which the catalog was generated.
                   Can be a string with the h5py key or an integer
                   with the snapshot number. If snapshot=None then the
                   file will be treated as an old IRATE.

        catalog - The AHF or Rockstar catalog. Can be a string with 
                  the h5py key or simply 'AHF'/'Rockstar', if given as
                  the later the first AHF/Rockstar catalog found will
                  be used   

        halo_id - The halo ID in the catalog, for Rockstar catalogs 
                  this alone identifies the halo. For AHF catalogs
                  IDs are not unique across MPI tasks and the 
                  halo_arg parameter may be necesary to identify
                  the desired halo

        halo_arg - The argument given to the catalog data arrays   
                   to identify the wanted halo, 
                   e.g. Vmax = Catalog['Vmax'][...][halo_arg]
  
        '''

        if (halo_id == None) and (halo_arg == None):
            errmsg = 'Either the halo_id or a halo_arg needs to be '\
                'given in order to initialize a new halo object\n'
            errmsg = errmsg + self.__init__.__doc__
            raise ValueError(errmsg)
                    
        # set some given attributes
        self.irate_file = irate_file
        if type(snapshot) == int: snapshot = 'Snapshot%.5d' % snapshot
        self.snapshot = snapshot
        self.id = halo_id
                        
        # set simulation properties
        self.sim_props = irate_file_mod.translate_filename_hyades(irate_file)
        
        # open snapshot
        irate = h5py.File(irate_file,'r')
        if snapshot:
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
            
        self.catalog = catalog
        self.catalog_attrs = {}
        # pass catalog attributes
        for key in C.attrs.keys():
            self.catalog_attrs[key] = C.attrs[key]

        # Find halo in catalog
        if halo_arg != None:
            arg = halo_arg
        else:    
            try:
                ID =  C['ID'][...]
                arg = argwhere(ID == halo_id)[:,0]
                if len(arg) > 1:
                    errmsg = "The halo ID's in this catalog are not unique try "\
                        "using halo_arg instead to identify halo\n"
                    errmsg = errmsg + self.__init__.__doc__
                    raise ValueError(errmsg)
                elif not arg:
                    errmsg = 'There is no halo with halo_id = %d in snapshot %s '\
                        'of catalog %s' % (halo_id, snapshot, catalog) 
                    raise ValueError(errmsg)
                else:
                    arg = arg[0]
            except KeyError:
                errmsg = 'No ID data found in catalog, most probably because ' \
                    'it is an old AHF catalog. use halo_arg instead to '\
                    'identify halo.\n'
                errmsg = errmsg + self.__init__.__doc__
                raise ValueError(errmsg)
        self.cat_arg = arg
        
        # Copy and set halo properties straight from catalog
        self.props = {}
        self.units = {}
        for key in C.keys():
            if ('RadialProfiles' not in key) and ('ParticleData' not in key):
                self.props[key] = C[key][...][arg]
                try: unitname = C[key].attrs['unitname']
                except KeyError: unitname = None
                try: unitcgs = C[key].attrs['unitcgs']
                except KeyError: unitcgs = None                
                self.units[key] = {'unitname': unitname, 'unitcgs': unitcgs}
                if key == 'ID':
                    self.id = C[key][...][arg]

        # Close the irate file
        irate.close()        


    def print_properties(self):
        '''Print some properties of this halo'''
        
        h = visnap.h
        print 'halo properties:'
        print 'Vmax = %g km/s, Rmax = %g kpc' % (self.props['Vmax'], 
                                                 self.props['Rmax']/h) 
        print 'Mvir = %g Msun, Rvir = %g kpc' % (self.props['Mvir']/h,
                                                 self.props['Rvir']/h)
        print 'Center [Mpc] = ', self.props['Center'] 
        print 'Velocity [km/s] = ', self.props['Velocity']
        print 'Npart = ', self.props['npart']
        if 'AHF' in self.catalog: 
            print 'fMhires = ', self.props['fMhires']
            print 'com_offset = ', self.props['com_offset'] 
            print 'mbp_offset = ', self.props['mbp_offset']


    def track(self, trees_path='./trees/', trees_file_base='tree', ncpus='all' ):
        '''
        Find this halo tree and trace its properties back in time. This only
        works for Rockstar halos.

        Input:
         
         trees_path - The path to the merger trees files

         trees_file_base - The base name of the trees files

         ncpus - If 'all' all the avilable cores in the node will be used, else set
                 to the number of cores that you want to use
         
        Output:

          halo_past_props - A dictionary containing the halo properties over time,
                            all halo properties in the rockstar tree are included 
                            with the same key names. 
        '''
        from visnap.general.halo_track import find_tree, halo_track
        # Check if this halo comes from a Rockstar catalog
        if not 'Rockstar' in self.catalog: 
            print 'At the moment halo tracking only works for halos found using '\
                'Rockstar'
            sys.exit()
        
        # Get needed properties
        halo_id = self.id    
        scale = round(self.catalog_attrs['a'], 4)
        
            
        # Find tree
        try:
            self.tree, self.depth_first_id = find_tree(halo_id, scale, trees_path,
                                                       trees_file_base, ncpus)
        except IOError:
            try:
                trees_file_base = self.irate_file.split('/')[-1]
                trees_file_base = trees_file_base.rstrip('irate.hdf5')+'trees'
                self.tree, self.depth_first_id = find_tree(halo_id, scale,
                                                           trees_path,
                                                           trees_file_base, ncpus)
            except IOError:
                print 'The merger trees files were not found, use the trees_path '\
                    'and trees_file_base inputs to give the correct location and '\
                    'base name of the trees files'
                sys.exit()
                
        
        # Track
        self.past_props = halo_track(self.tree, self.depth_first_id,
                                     trees_path, trees_file_base)
        return self.past_props
            

    def get_subhalos(self,  vcut_factor=0.05, rcut_factor=1.0, min_npart=100):
        '''
        Find the subhalos of this halo and create list of halo objects from
        them, make this list available for later use as self.subhalos 

        Input

         vcut_factor - The lowest subhalo vmax to consider in units of the 
                       host Vmax

         rcut_factor - Only subhalos with r < rcut_factor*host_Rvir will be
                       returned
   
         min_npart - Only subhalos with at least min_npart number of particles
                     will be included 
                  
        Output

         subhalo_list - A list of halo objects with all the found subhalos
        '''
        from visnap.general.find_halos import find_subhalos
        
        subhalo_list = find_subhalos(self, vcut_factor, rcut_factor,
                                     min_npart) 
        self.subhalos = subhalo_list

        return subhalo_list    
       
     
    def get_profiles(self, particles_file=None, Nbins=50, Ncpus='all',
                     remove_subs=False, project_along=None, only2d=False): 
        '''
        Get the radial profiles of this halo and set them up
        as attributes for future use

        Input:

         particles_file - The rockstar particles files containing the particle
                          information  for the halo with ID = rockstar_halo_id.
                          If None it will be searched based on the IRATE filename
     
         Nbins - Nuber of bins to use when calculating profiles
                 for Rockstar halos

         Ncpus - For rockstar catalogs the particle data can be
                 extracted in parallel. Set ncpus to 'all'
                 to use all the available cores in the node, else
                 set to the number of cores you want to use
         
         remove_subs - If True the particles of the subhalos will be removed
                       and only the smooth halo distribution will be used (Only
                       works for Rockstar catalogs at the moment)

         project_along - If you want 2d projection profiles you can set
                         project_along to a vector defining the projection axis,
                         e.g. project_along=(1,0,0) will calculate 2d
                         projection profiles along the x-axis (Only works for
                         Rockstar catalogs at the moment)

         only2d - If set to True only 2d projection profiles will be produced  
              
        '''
        from visnap.functions.mis import find_rpower

        # open snapshot
        irate = h5py.File(self.irate_file,'r')
        snap = irate[self.snapshot]

        # open catalog
        C = snap[self.catalog]
        
        # get profiles
        profiles = {}
        if 'AHF' in self.catalog:
            P = C['RadialProfiles']
            rhob = C.attrs['rho_back(z)']/10**9 # [Msun/kpc^3 h^2]
            for key in P.keys():
                profiles[key] = P[key][...][self.cat_arg]
                if (key == 'dens') or (key == 'ovdens'):
                    profiles[key] *= rhob
            if project_along!=None:
                print '2d projection profiles are only supported for Rockstar '\
                    'catalogs at the moment'
        else:
            import visnap.general.halo_particles as halo_particles
            print 'Attempting to calculate profiles from Rockstar '\
                'particle data, this could take some time when '\
                'done for the first time'
            
            pfile_ext = str(int(self.snapshot.strip('/Snapshot'))) + '.particles'
            particles_file = self.irate_file.replace('irate.hdf5',pfile_ext)
            pdata = halo_particles.get_rockstar_halo_particles(self.id, particles_file,
                                                   snap, Ncpus)
            if remove_subs:
                asign_id, int_id =  pdata[:,7], pdata[:,8]
                pdata = pdata[asign_id==int_id] 
                print 'Removing particles belonging to subhalos, only %d ' \
                    'particles forming the smooth halo distribution will be used ' \
                    '%(pdata[:,0].size)'
            
            halo_center = self.props['Center']
            halo_velocity = self.props['Velocity']
            if 'Mpc' in C['Center'].attrs['unitname']:
                #convert particle positions Mpc/h (comoving) -> Kpc/h (physical)
                pdata[:,0:3] = pdata[:,0:3].astype(float)*1000.0*C.attrs['a'] 
                #convert halo position Mpc/h (comoving) -> Kpc/h (physical)
                halo_center = halo_center*1000.0*C.attrs['a']
           
            mpdm = self.sim_props['mpdm']
            
            if not only2d:    
                print 'Calculating profiles from particle data'    
                cp = halo_particles.calculate_profiles(pdata[:,0:6], halo_center,
                                           halo_velocity, Nbins)
                r_p, Nshell_p, Nenclosed_p, dens_p, avgDens_p, Vdisp_p = cp
                profiles = {'r': r_p, 'npart_shell': Nshell_p, 'npart': Nenclosed_p,
                            'M_in_r': Nenclosed_p*mpdm, 'dens': dens_p*mpdm, 
                            'ovdens': avgDens_p*mpdm, 'sigv': Vdisp_p}
                self.profiles = profiles
            
            if project_along!=None:
                print 'Calculating 2d profiles from particle data'
                print 'Projection axis = ', project_along
                cp = halo_particles.calculate_2dprofiles(pdata[:,0:6], halo_center,
                                            halo_velocity, project_along, Nbins)
                r_2dp, Nshell_2dp, Nenclosed_2dp, dens_2dp, avgDens_2dp,\
                    Vlosdisp_2dp, Vtandisp_2dp = cp
                profiles2d = {'R': r_2dp, 'npart_ring': Nshell_2dp, 'npart': Nenclosed_2dp,
                              'M_in_R': Nenclosed_2dp*mpdm, 'Sdens': dens_2dp*mpdm, 
                              'avgSdens': avgDens_2dp*mpdm, 'sigVlos':  Vlosdisp_2dp,
                              'sigVtan': Vtandisp_2dp}
                self.profiles2d = profiles2d
                self.profiles2d_projection_axis = project_along
            
        # set calculated profiles properties
        if not only2d:        
            r, npart,= abs(profiles['r']), profiles['npart'] 
            ovdens, dens = profiles['ovdens'], profiles['dens']
            #Resolved radii based on Power et al. Relaxation time scale (Eq. 20) 
            rPower, argRes, argNotRes = find_rpower(r, npart, ovdens)
            self.profiles_props = {'rPower': rPower, 'argRes': argRes, 
                                   'argNotRes': argNotRes}
            if self.sim_props['dmName'] == 'SIDM':
                mpdm, hsi, epsilon = [self.sim_props['mpdm'], self.sim_props['hsi'],
                                      self.sim_props['epsilon']]
                arg = argwhere(dens[argRes] > mpdm*(0.2/(hsi*2.8*epsilon))**3)
                if len(arg) > 0:
                    self.profiles_props['rmax_sidm'] = r[argRes][arg].max()
                else:
                    self.profiles_props['rmax_sidm'] = 0.0
            
        # Close the irate file
        irate.close()

    def plot_density(self):
        '''Plot density profile'''
        import visnap.plot.halo_profiles as halo_profiles
        halo_profiles.density_profile(self)


        
    

        
        
