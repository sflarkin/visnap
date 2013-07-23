'''
VISNAP (irVIne coSmology aNAlysis Package) is a package of python
modules for the analysis of simulation data saved under the IRATE (IRvine 
Astorphysical simulaTion structurE) format.

VISNAP is being developed by the UC Irvine cosmology group, for questions
contact: Miguel Rocha - rocham@uci.edu, Shea Garrison-Kimmel - sgarriso@uci.edu
or Jose Onorbe - jonorbeb@uci.edu.
'''

# Define belwow any parameters that are to be used package-wide,
# they are available to all VISNAP modules by calling visnap.parameter

#UNITS

UnitLength_in_cm = 3.085678e21  # kpc
UnitMass_in_g = 1.989e33 # Msun

#PHYSICS AND COSMOLOGY (WMAP7)

G = 4.30071e-6 # [Msun^-1 kpc (km/s)^2]
omega_matter = 0.266
rho_b = 73.83  #background density at z=0 [Msun/kpc^3 h^2]
rho_crit = rho_b/omega_matter # at z=0 
delta_ovdens = 365.7 # at z=0
h = 0.71
dh = 1.9e-5
keq = 0.0732*omega_matter*h**2

#FLAGS


#EVERYTHING ELSE


 
