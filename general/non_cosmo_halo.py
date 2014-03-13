'''Module containing tools for generating halo objects from idealized
(non-cosmo) runs'''

import sys
import h5py
from numpy import *
import visnap
#import pdb

class new_halo:
    '''Create a new halo object'''

    def __init__(self, irate_file, snapshot, halo_id=None, center='pot'):
        '''
        Initialize new halo object with basic attributes and properties

        Input:
        
         irate_file - The name of the IRATE file containing the catalog

         snapshot - The snapshot from where to extract the DM particles.
                    Can be a string with the h5py key or an integer
                    with the snapshot number
                    
         halo_id - The ID given to this halo, Can be an integer or a string. If
                   None it will set as "<snapshot_name>_halo1"

         center - Either 'mean' to use the mean of the particle positions as the
                  the center of the halo, or 'pot' for the minimum of the potential.
                  Also you can set it as (x, y, z) or [x, y, z].
        '''           

        # set some given attributes
        self.irate_file = irate_file
        if type(snapshot) == int: snapshot = 'Snapshot%.5d' % snapshot
        self.snapshot = snapshot
        if halo_id: self.id = halo_id
        else: self.id = snapshot+'\_halo1'      
        
        # get profiles and profile derived quantities
        self.get_profiles(center=center)

    def print_properties(self):
        '''Print some properties of this halo'''
        
        print 'halo properties:'
        print 'Vmax = %g km/s, Rmax = %g kpc' % (self.props['Vmax'], 
                                                 self.props['Rmax']) 
        print 'Center [Mpc] = ', self.props['Center'] 
        print 'Velocity [km/s] = ', self.props['Velocity']
        print 'Npart = ', self.props['npart']
        

    def get_profiles(self, Nbins=50, center='pot'):  
        '''
        Get the radial profiles of this halo and set them up
        as attributes for future use

        Input:
      
         Nbins - Number of bins to use when calculating profiles

         center - Either 'mean' to use the mean of the particle positions as the
                  the center of the halo, or 'pot' for the minimum of the potential.
                  Also you can set it as (x, y, z) or [x, y, z].
                 
        '''

        from visnap.functions.mis import find_rpower
        import visnap.general.halo_particles as hp

        # Open snapshot
        irate = h5py.File(self.irate_file,'r')
        snap = irate[self.snapshot]
        Dark_Halo = snap['ParticleData']['Dark_Halo']
        
        # Get particle data
        Mass = Dark_Halo['Mass'][...]
        Pos = Dark_Halo['Position'][...]
        try:
            Pot = Dark_Halo['Potential'][...]
        except KeyError:
            center = 'mean'
            print 'Warning: There is no field for the potential of the '\
                'particles, the center of the halo will be set as the mean '\
                'of the particle positions.'
            
        Vel = Dark_Halo['Velocity'][...]
               
        # Find center and bulk velocity
        if center=='pot':
            halo_center = Pos[argwhere(Pot == Pot.min())[0,0]]
        elif center=='mean':
            halo_center = array([mean(Pos[:,0]), mean(Pos[:,1]), mean(Pos[:,2])])
        else:
            halo_center = array(center)
            
        halo_velocity = array([Vel[:,0].mean(), Vel[:,1].mean(),
                               Vel[:,2].mean()]) 
        self.props = {'Center': halo_center, 'Velocity': halo_velocity,
                      'npart': len(Mass)}

        # Calculate profiles
        if len(set(Mass)) != 1:
            print 'There are particles with different masses, for now this is'\
                ' not supported'
            sys.exit()
        pdata = array([Pos[:,0], Pos[:,1], Pos[:,2], Vel[:,0], Vel[:,1],
                       Vel[:,2]])
        pdata = pdata.transpose()
        cp = hp.calculate_profiles(pdata, halo_center, halo_velocity, Nbins)
        r_p, Nshell_p, Nenclosed_p, dens_p, avgDens_p, Vdisp_p = cp
        mpdm = Mass[0]*1e10
        print 'Warning: The mass of particles is assumed to be in Gadget units'\
            ' (i.e. 1e10 Msun), you will need to correct if different.'
        
        profiles = {'r': r_p, 'npart_shell': Nshell_p, 'npart': Nenclosed_p,
                    'M_in_r': Nenclosed_p*mpdm, 'dens': dens_p*mpdm, 
                    'ovdens': avgDens_p*mpdm, 'sigv': Vdisp_p}
        
        # Set profiles and calculated profiles properties
        self.profiles = profiles
        r, npart = abs(profiles['r']), profiles['npart'] 
        ovdens, dens = profiles['ovdens'], profiles['dens']
        # Resolved radii based on Power et al. Relaxation time scale (Eq. 20) 
        rPower, argRes, argNotRes = find_rpower(r, npart, ovdens)
        self.profiles_props = {'rPower': rPower, 'argRes': argRes, 
                              'argNotRes': argNotRes} 
        
        # Find Vmax and Rmax
        G = visnap.G
        Menc = profiles['M_in_r']
        vcirc = sqrt(G*Menc/r)
        self.props['Vmax'] = max(vcirc)
        self.props['Rmax'] = r[argwhere(vcirc == max(vcirc))]
       
        self.sim_props = {'mpdm': mpdm, 'dmName': 'nonCosmo'}

        # Close the irate file
        irate.close()


    def plot_density(self):
        '''Plot density profile'''
        import visnap.plot.halo_profiles as hp
        hp.density_profile(self)
 
