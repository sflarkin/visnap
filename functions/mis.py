'''Module containing miscellaneous short functions'''

from numpy import *
import visnap

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


def straight_line(x, x_cross, y_cross, dx=None, dy=None, slope=1):
    '''
    For the given x return the corresponding y in the straight line
    defined by x_cross, y_cross and slope (slope defaults to 1). 
    If dx and dy are given then slope = dy/dx.
    '''
    if (dx!=None) & (dy!=None): slope = dy/dx
    a = y_cross-slope*x_cross
    return a + slope*x

