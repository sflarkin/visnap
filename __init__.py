'''
VISNAP (irVIne coSmology aNAlysis Package) is a package of python
modules for the analysis of simulation data saved under the IRATE (IRvine 
Astorphysical simulaTion structurE) format.

VISNAP is being developed by the UC Irvine cosmology group, for questions
contact: Miguel Rocha - rocham@uci.edu, Shea Garrison-Kimmel - sgarriso@uci.edu
or Jose Onorbe - jonorbeb@uci.edu.
'''
import numpy as np

# Define belwow any parameters that are to be used package-wide,
# they are available to all VISNAP modules by calling visnap.parameter

#UNITS

UnitLength_in_cm = 3.085678e21  # kpc
UnitMass_in_g = 1.989e33 # Msun

#PHYSICS

G = 4.30071e-6 # [Msun^-1 kpc (km/s)^2]
c = 299792.458 # [km/s]

#COSMOLOGY (WMAP7+BAO+H0, Komatsu et al. 2011) 

h = 0.702
omega_matter = 0.272
omega_lambda = 0.728
omega_baryons = 0.0455
ns = 0.961
sigma_8 = 0.807
rho_crit = 3.0*(100.0/1000.0)**2/(8.0*np.pi*G) # at z=0 [Msun/kpc^3 h^2]
rho_b = rho_crit*omega_matter  # background density at z=0  [Msun/kpc^3 h^2]
delta_ovdens = 360.6 # at z=0
dh = 1.9e-5
keq = 0.0732*omega_matter*h**2

#FLAGS


#EVERYTHING ELSE


 
