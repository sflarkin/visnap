'''Module containing model functions for radial profiles'''

from numpy import *

def rhoNFW(logr, param):
    '''NFW density profile,  param = (normalization, scale_radius)'''
    r = 10**logr
    rhos,rs = param[0],param[1]
    return log10(rhos/((r/rs)*(1+r/rs)**2))

def MNFW(logr, param):
    '''NFW mass profile,  param = (normalization, scale_radius)'''
    r = 10**logr
    rhos,rs = param[0],param[1]
    return log10(4.0*pi*rhos*rs*rs*rs*(-r/(r+rs) + log((r+rs)/rs)))

def dMNFWdrs(r, param):
    '''NFW dM/d_rs derivative,  param = (normalization, scale_radius)'''
    rhos,rs = param[0],param[1]
    fact1 = 4.0*pi*rhos*rs*rs
    fact2 = -r*(4.0*r + 3.0*rs)/((r+rs)*(r+rs))
    fact3 = 3.0*log((r+rs)/rs)
    return fact1*(fact2+fact3)

def dMNFWdrhos(r, param):
    '''NFW dM/d_rhos derivative,  param = (normalization, scale_radius)'''
    rhos,rs = param[0],param[1]
    fact1 = 4.0*pi*rs*rs*rs
    fact2 = -r/(r+rs)
    fact3 = log((r+rs)/rs)
    return fact1*(fact2+fact3)

def SurfaceDensityNFW(logR, param):
    '''NFW surface density profile, param = (1d density normalization, scale_radius)'''
    R = 10**logR
    rhos,rs = param[0],param[1]
    x = R/rs
    sigmas = 2*rs*rhos
    if x < 1:
        fac1 = sigmas/(x*x - 1.0)
        fac2 = 1.0-2.0*arctanh(sqrt((1.0-x)/(1.0+x)))/sqrt(1-x*x)
        sigma = fac1*fac2
    elif x == 1:
        sigma = sigmas/3.0
    else:
        fac1 = sigmas/(x*x - 1.0)
        fac2 = 1.0-2.0*arctanh(sqrt((x-1.0)/(1.0+x)))/sqrt(x*x-1)
        sigma = fac1*fac2
    return log10(sigma)

def AvgSurfaceDensityNFW(logR, param):
    '''NFW mean enclosed surface density profile, param = (1d density normalization, scale_radius)'''
    R = 10**logR
    rhos,rs = param[0],param[1]
    x = R/rs
    sigmas = 4*rs*rhos
    if x < 1:
        fac1 = sigmas/(x*x) 
        fac2 = 2.0*arctanh(sqrt((1.0-x)/(1.0+x)))/sqrt(1-x*x)
        fac3 = log(x/2.0)
        sigma = fac1*(fac2+fac3)
    elif x == 1:
        sigma = sigmas*(1+log(0.5))
    else:
        fac1 = sigmas/(x*x)
        fac2 = 2.0*arctanh(sqrt((x-1.0)/(1.0+x)))/sqrt(x*x-1)
        fac3 = log(x/2.0)
        sigma = fac1*(fac2+fac3)
        
    return log10(sigma)    

def rhoBurkert(logr, param):
    '''Burkert density profile,  param = (normalization, core_radius)'''
    r = 10**logr
    rhos,rc = param[0],param[1]
    return log10(rhos/((1+(r/rc)**2)*(1+r/rc)))  

def rhoModBurkert(logr, param):
    '''Burkert density profile with 2 scale radii, param = (normalization,
    core_radius, scale_radius'''
    r = 10**logr
    rhos, rc, rs = param[0], param[1], param[2]
    return log10(rhos/((1+(r/rs)**2)*(1+r/rc)))  

def rhoZavala(logr, param):
    '''Another Burkert like density profile with 2 scale radii, param = (normalization,
    core_radius, scale_radius'''
    r = 10**logr
    rhos, rc, rs = param[0], param[1], param[2]
    return log10(rhos*rs**3/((r+rc)*(r+rs)**2))  
                 
def MBurkert(logr, param):
    '''Burkert mass profile,  param = (normalization, core_radius)'''
    r = 10**logr
    rhos,rc = param[0],param[1]
    return log10(pi*rhos*rc*rc*rc*(-2.0*arctan(r/rc) + log((r+rc)*(r+rc)*(r*r+rc*rc)/rc**4)))

def dMBurkertdrc(r, param):
    '''Burkert dM/d_rc derivative,  param = (normalization, core_radius)'''
    rhos,rc = param[0],param[1]
    fact1 = pi*rhos*rc*rc
    fact2 = -4.0*r*r*r
    fact3 = 3.0*(r+rc)*(r*r+rc*rc)*(-2.0*arctan(r/rc)+log((r+rc)*(r+rc)*(r*r+rc*rc)/rc**4))
    fact4 = (r+rc)*(r*r+rc*rc)
    return fact1*(fact2+fact3)/fact4

def dMBurkertdrhos(r, param):
    '''Burkert dM/d_rhos derivative,  param = (normalization, core_radius)'''
    rhos,rc = param[0],param[1]
    fact1 = pi*rc*rc*rc
    fact2 = -2.0*arctan(r/rc)
    fact3 = log((r+rc)*(r+rc)*(r*r+rc*rc)/rc**4)
    return fact1*(fact2+fact3)          

def rhoPsudoIso(logr, param):
    '''Pseudo-Isothermal density profile, param = (normalization, core_radius)'''
    r = 10**logr
    rhos,rc = param[0],param[1]
    return log10(rhos/(1+(r/rc)**2)) 

def MPsudoIso(logr, param):
    '''Pseudo-Isothermal mass profile, param = (normalization, core_radius)'''
    r = 10**logr
    rhos,rc = param[0],param[1]
    return log10(4.0*pi*rhos*rc*rc*rc*(r/rc-arctan(r/rc)))

def dMPsudoIsodrc(r, param):
    '''Pseudo-Isothermal dM/d_rc derivative, param = (normalization, core_radius)'''
    rhos,rc = param[0],param[1]
    fact1 = 4.0*pi*rhos*rc
    fac2 = 3.0*r - r*r*r/(r*r+rc*rc)
    fact3 = -3.0*rc*arctan(r/rc)
    return fact1*(fact1+fact2)

def dMPsudoIsodrhos(r, param):
    '''Pseudo-Isothermal dM/d_rhos derivative, param = (normalization, core_radius)'''
    rhos,rc = param[0],param[1]
    fact1 = 4.0*pi*rc*rc
    fact2 = r-rc*arctan(r/rc)
    return fact1*fact2
        
def rhoEinasto(logr, param):
    '''Einasto density profile, param = (normalization, scale_radius, alpha )'''
    r = 10**logr
    rhos,rs,alpha = param[0],param[1],param[2]
    return log10(rhos*exp(-(r/rs)**alpha))

def rhoModNFW(logr, param):
    '''Modified NFW density profile, param = (normalization, scale_radius, inner_slope)'''
    r = 10**logr
    rhos,rs,c = param[0],param[1],param[2]
    return log10(rhos/((r/rs)**c*(1+r/rs)**(3-c)))

def rhoManojs(logr, param):
    '''Suggested density profile for SIDM by Manoj's,  param = (normalization, scale_radius, core_radius)'''
    r = 10**logr
    rhos,rs,rc = param[0],param[1],param[2]
    return log10(rhos/((r/rs)*(1+r/rs)**2 + rc/rs))


