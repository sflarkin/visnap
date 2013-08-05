'''
CentroidInfer.py made from template.py
Author: Miguel Rocha <mrocha@gmail.com>
Mostly put together by: Karen Y. Ng <karenyng@ucdavis.edu>
Purposes: 
gets the projected positions and radial velocities of the subhalos (galaxies)
then calls Will 's existing modules for 
-creating a fits map 
-creating bootstrap analysis of the projected centroid based on biweight stats.
-estimating galaxy centroids

prequisites by William A. Dawson:
CAT.py 
CenterVarianceBoot.py 
input: 
irate_files 

output:
- a galaxy density fits map   

License: BSD
'''
from __future__ import division

#! /usr/bin/python   
import numpy as np
import math
import os, sys, time
import matplotlib.pyplot as plt
from glob import glob
from optparse import OptionParser
import pdb

def ang_coord(subhalo1,obs_pos):
    '''
    Purposes: computes the line of angular coordinates based on
    position of the observer 
    input: 
    subhalo1 = subhalo object returned by zoom_halo.get_subhalo()
    obs_pos = numpy array that contains the x,y,z position of the observer
    output: 
    RA = RA in degrees 
    DEC = DEC in degrees 
    '''
    subh_center = subhalo1.props['Center']
    subh_center_sigma=subhalo1.props['PositionUncertainty']
    v3d_sigma = subhalo1.props['VelocityUncertainty']
    #have to think about if we want RA increase to +x or -x 
    del_center = subh_center-obs_pos
    del_r = np.linalg.norm(del_center)
    delta = 180./np.pi*np.arcsin(del_center[1]/del_r) 
    #have z=0 to be where RA = 0 
    alpha = 180./np.pi*np.arctan(del_center[0]/del_center[2]) 

    return alpha,delta 

def proj_v_rad(subhalo1,obs_pos):
    '''
    This function returns the projected radial velocity of a particular subhalo
    relative to the observer, that is due to peculiar motion
    inputs: 
    subhalo = subhalo object returned by halo.find_subhalos
    obs_pos = numpy array with x,y,z position of the observer
    output:
    v_rad = radial velocity with +ve value moving away from us 
    and -ve value moving towards us 
    '''
    v3d = subhalo1.props['Velocity']
    subh_center = subhalo1.props['Center']
    #have to think about if we want RA increase to +x or -x 

    del_center = subh_center-obs_pos
    #compute radial velocities by first computing the unit vector:
    subh_ucenter = (subh_center - obs_pos)/\
        np.sqrt(np.dot(subh_center-obs_pos,subh_center-obs_pos))
    # project the 3d velocity unto the radial direction
    v_rad = np.dot(v3d,subh_ucenter)*subh_ucenter
    #divide by the unit vector in order to get the sign right
    v_rad  = (v_rad / subh_ucenter)[0]                                        
    if np.isnan(v_rad)==True:
        print 'visnap.centroidInfer.proj_v_rad:\n'
        print 'NAN value detected'
        sys.exit()

    return v_rad 

if __name__ != '__main__':
    print '%s is only script, there is nothing to import from it.' % __name__
    sys.exit()

#Parsing Options

parser = OptionParser()

parser.add_option("-f", "--files", dest="irate_files",
                  default = "",
                  help="IRATE files to read, e.g. <path/filename> for a single " \
                  "file and <path1/file1,path2/file2,...> or <[path1/file1, "\
                  "path2/file2, ...]>  for multiple files. Leave blank to read "\
                  "all files with BASE*.EXT on PATH, where BASE,EXT and PATH " \
                  "are defined with the -b,-e and -p options respectively")
parser.add_option("-b", "--base", dest="base",
                  default="", 
                  help="Base of IRATE files (default = '')")
parser.add_option("-e", "--extension", dest="ext",
                  default="_irate.hdf5", 
                  help="Extension of IRATE files (default = _irate.hdf5)")
parser.add_option("-p", "--path", dest="path",
                  default="catalogs/", 
                  help="Path to IRATE files (default = ./catalogs). NOTE: if "\
                  "files are given with the -f option it is assumed that the "\
                  "paths are included with the names, thus the -p parameter "\
                  "will be ignored.")

(options, args) = parser.parse_args()

irate_files = options.irate_files
base = options.base
ext = options.ext
path = options.path

#Set the irate_files variable
if irate_files == "":
    if '+' in ext:
        ext = ext.split("+")
        irate_files = [glob(path + base + '*' + thisext) for thisext in ext]
        irate_files = concatenate(irate_files)
    else:    
        irate_files = glob(path + base + '*' + ext)
    irate_files.sort()
elif "[" in irate_files:
    irate_files = irate_files.replace("[","")
    irate_files = irate_files.replace("]","")
    irate_files = irate_files.replace(" ","")
    irate_files = irate_files.split(",")
    irate_files.sort()  
elif "," in irate_files:
    irate_files = irate_files.replace(" ","")
    irate_files = irate_files.split(",")
    irate_files.sort()
else:
    irate_files = [irate_files,]

print "\nThe following %d IRATE files were parased for analysis:\n"%len(irate_files),irate_files[:]
print '' 

import visnap.general.translate_filename as translate_filename
import visnap.general.find_halos as find_halos
import visnap.plot.halo_profiles as plot_profiles

#-----initialize some variables--------------------
colnames = ('ra','delta','v_rad','v_rad_sigma','centerx','centery','centerz',
            'center_sigma','vmax')
write_out=True
prefix = 'gal_cat_xz-plane'
c = 3.e5 #in km/s
rabin = 200
N_boot = 100
out_path = '/home/karen/ResearchCode/centroid/'

#-----start code --------------------
halos = []
for irate_file in irate_files:
    #simulation properties 
    simprops = translate_filename.hyades(irate_file)
    zoom_id = simprops['zoom_id']
    
    #can find the snapshot number by doing $h5ls irate.hdf5 filename
    zoom_halo, dist = find_halos.find_zoom_halo('../650Box_clusters_zooms.txt',
                                                zoom_id,irate_file,
                                                'Snapshot00148','Rockstar')
    halos.append(zoom_halo)
    
    halo1 = halos[0] 
    #try to put the location of the observer at different projected views 
    obs_posx = halo1.props['Center'][0]#0.
    obs_posy = 0. #halo1.props['Center'][1]#0.
    obs_posz = halo1.props['Center'][2]#0.
    obs_pos = [obs_posx,obs_posy,obs_posz]
    h0 = halo1.catalog_attrs['h0']
    Ol = halo1.catalog_attrs['Ol']
    Om = halo1.catalog_attrs['Om']
    a = halo1.catalog_attrs['a']
    haloRA,haloDEC = ang_coord(halo1,obs_pos) 
    halo_v_rad = proj_v_rad(halo1,obs_pos)

    subhalo_list = halo1.get_subhalos(1.5)
    if write_out==True:
        iprefix1,iprefix2,junk=irate_file.split('.')
        catalog = prefix+'.txt' 
       
        print 'opening file: '+catalog+' for catalog output'
        f = open(out_path+catalog,'w')
        f.write('#visnap.scripts.centroidInfer.py output:\n')
        f.write('#This catalog is written using inputs from:\n')
        f.write('#'+irate_file+'\n')

        f.write('#The host halo properties are \n'+\
        '#RA Dec v_sigma centerx centery centerz center_sigma vmax rvir mvir\n'+\
        '#\t{0}\t {1}\t {2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10:1.5e}\n'.format(haloRA,
                                        haloDEC,
                                        halo_v_rad,
                                        halo1.props['VelocityUncertainty'],
                                        halo1.props['Center'][0],
                                        halo1.props['Center'][1],
                                        halo1.props['Center'][2],
                                        halo1.props['PositionUncertainty'],
                                        halo1.props['Vmax'],
                                        halo1.props['Rvir'],
                                        halo1.props['Mvir']))
                                             
        f.write('#based on the observer at :{0}\n'.format(obs_pos))
        for i in range(len(colnames)):
            f.write('#ttype'+str(i)+' = '+colnames[i]+'\n')
    
    for i in range(len(subhalo_list)):
        subhalo1= halo1.subhalos[i]
        alpha, delta = ang_coord(subhalo1, obs_pos)
        v_rad = proj_v_rad(subhalo1, obs_pos) 
        #ignore projection effects from position uncertainty
        v_rad_sigma = subhalo1.props['VelocityUncertainty']
        #----writing out the outputs -----
        if write_out ==True:
            f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(alpha,
                                                      delta,
                                                      v_rad,
                                                      v_rad_sigma,
                                                      subhalo1.props['Center'][0],
                                                      subhalo1.props['Center'][1],
                                                      subhalo1.props['Center'][2],
                                                      subhalo1.props['PositionUncertainty'],
                                                      subhalo1.props['Vmax']))
    if write_out==True:
        f.close()

