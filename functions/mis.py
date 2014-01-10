'''Module containing miscellaneous short functions'''

from numpy import *
import visnap


def particle_mass(box_size, res_level):
    '''
    Find the particle mass in a simulation

    Input:

     box_size - 1d side lenght of simulation box in [Mpc/h] 

     res_level - Resolution level, i.e 2^(res_level) = Nparticles^(1/3)
     
    Output:
    
      mdark, mgas - Mass of dark matter and gas particles in [Msun/h]
    '''
    
    rho_b = visnap.rho_b
    Ob = visnap.omega_baryons
    mp = (rho_b*(box_size*1000.0)**3)/2.0**(3.0*res_level)
    mgas = mp*Ob

    mdark = mp-mgas
    return mdark, mgas


def optimal_resolution(box_size, res_level, Mvir, Rvir, dark_matter_only=True):
    '''
    Find the optimal resolution to simulate a given dark matter halo following
    power et al. 2003 criterion, i.e. opt_res = 0.5*4*R200/sqrt(N200)~1.6 Rvir/sqrt(Nvir)

    Input:
    
     box_size - 1d side lenght of simulation box in [Mpc/h]

     res_level - Resolution level, i.e 2^(res_level) = Nparticles^(1/3)

     Mvir - Virial mass of dark matter halo in [Msun/h]

     Rvir - Virial radius of dark matter halo

     dark_matter_only - Set to false if the simulation includes baryons

    Output:
    
     opt_res - optimal resolution in [same units as Rvir]   
    '''

    mdark, mgas = particle_mass(box_size, res_level)
    if dark_matter_only:
        mp = mdark+mgas
    else:
        mp = mdark
  
    opt_res = 1.6*(Rvir/sqrt(Mvir/mp))
    return opt_res

    
def find_rpower(r, nenc, avg_dens):
    '''
    Find the radius at which the numerical two body scattering time-scale is
    comparable to the Hubble time.  Eq. 20 in power et al. 2003 is used.

    Input:

     r - radial bins

     nenc - number of enclosed particles at r

     avg_dens - average enclosed density at r

    Output:

     r_power - minimum resolved radius

     argRes - r[argRes] gives all r > r_power  

     argNotRes - r[argNotRes] gives all r < r_power
    '''
    
    rhoCrit = visnap.rho_crit
    argNotRes = argwhere(sqrt(200.0)*asfarray(nenc)*sqrt(rhoCrit/avg_dens)/(8.0*log(asfarray(nenc))) < 1.0)[:,0]
    r_power = r[argNotRes].max()
    argRes = argwhere(r > r_power)[:,0]
    return r_power, argRes, argNotRes


def findEinsteinR_NFW(rhos, rs, Zs, Zd):
    '''
    Find the Einstein radius for the given NFW halo

    Input:
    
     rhos - NFW Normalization in [Msun/Kpc^3]

     rs - Nfw scale radius in [Kpc]

     Zs - Redshift to lensing source 

     Zd - Redshift to deflector (i.e lens) 
          
    Output:

     Einstein radius in [Kpc]
    '''
    from scipy.optimize import brentq, fsolve
    from visnap.functions.radial_profiles import AvgSurfaceDensityNFW

    param = (rhos, rs)
    SigmaCrit = crit_surf_dens(Zs, Zd) 
    print 'Finding Einstain radius for NFW halo with rhos = %g and rs = %g' % param
    print 'Critical Surface Density [Msun/Kpc^2] = ', SigmaCrit
    print 'Average Surface Density at rs [Msun/Kpc^2] = ', 10**AvgSurfaceDensityNFW(log10(rs), param)
    sol = brentq(lambda R: AvgSurfaceDensityNFW(log10(R), param) - log10(SigmaCrit), 0.00001, 100*rs)                   
    return sol                

def crit_surf_dens(Zs, Zd):
    '''
    Find critical surface mass density 

    Input:
    
     Zs - Redshift to lensing source 

     Zd - Redshift to deflector (i.e lens) 
          
    Output:
    
     Sigma_crit - Critical surface density in [Msun/Kpc^2]
    '''

    G = visnap.G
    c = visnap.c
    Ds = AngDiamDist(Zs)*1000 # Mpc -> Kpc
    Dd = AngDiamDist(Zd)*1000
    Dds = AngDiamDist(Zs,Zd)*1000
    SigmaCrit = c*c*Ds/(4.0*pi*G*Dd*Dds)
    
    return SigmaCrit

def AngDiamDist(Z,Zo=0):
    '''
    Find the angular diameter distance for the given redshift z

    Input:
    
     Z - Redshift of observed object 
     
     Zo - Origin redshift (defaults to 0, i.e today)

    Output:

      D - Angular diameter distance in [Mpc]                
    '''
    from scipy.integrate import quad                   
    Om = visnap.omega_matter
    Ol = visnap.omega_lambda
    Ho = visnap.h*100.0
    c  = visnap.c                   
    a = 1.0/(1.0+Z)
    ao = 1.0/(1.0+Zo)
    
    D = quad(lambda ap: 1/(ap**2*Ho*sqrt(Om/ap**3 + Ol)), a, ao) 

    return c*a*D[0]      
                       

def straight_line(x, x_cross, y_cross, dx=None, dy=None, slope=1):
    '''
    For the given x return the corresponding y in the straight line
    defined by x_cross, y_cross and slope (slope defaults to 1). 
    If dx and dy are given then slope = dy/dx.
    '''

    if (dx!=None) & (dy!=None): slope = dy/dx
    a = y_cross-slope*x_cross
    return a + slope*x                       



