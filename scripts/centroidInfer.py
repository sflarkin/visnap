'''
CentroidInfer.py made from template.py
purposes: 
gets the projected positions and velocities of the subhalos (galaxies)
then calls Will 's existing modules for 
-creating bootstrap analysis
-estimating galaxy centroids

input: 
irate_files 

output:
- a galaxy density map in physical units from catalog  
- a galaxy density map that is returned by Will 's module

'''
from __future__ import division

#! /usr/bin/python   
import numpy as np
import math
import os, sys, time
from glob import glob
from optparse import OptionParser
#import pdb

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

print "\nThe following %d IRATE files were parsed for analysis:\n"%len(irate_files),irate_files[:]
print '' 

import visnap.general.translate_filename as translate_filename
import visnap.general.find_halos as find_halos
import visnap.plot.halo_profiles as plot_profiles

#-----initialize some variables--------------------
arg_list = ('objid','delta','rad','radv_sigma')
halos = []
write_out=True



#-----start code --------------------
for irate_file in irate_files:
    #simulation properties 
    simprops = translate_filename.hyades(irate_file)
    zoom_id = simprops['zoom_id']
    zoom_halo, dist = find_halos.find_zoom_halo('../650Box_clusters_zooms.txt',
                                                zoom_id, irate_file,
                                                'Snapshot00148', 'Rockstar')
    halos.append(zoom_halo)
    
#fig, ax, lines, legends = plot_profiles.density_profiles(halos
#simprop = visnap.general.translate_filename.hyades('GID_0650_12_1002_001_s01_h0.25_rockstar_irate.hdf5')

####only get subhalos within .3 rvir and with 1000 particles 
### to be changed later 
    halo1 = halos[0] 
    subhalo_list = halo1.get_subhalos(.3,1000)
    if write_out==True:
        iprefix1,iprefix2,junk=irate_file.split('.')
        f = open('/home/karen/ResearchCode/centroid/gal_cat_'+iprefix1+\
             '_'+iprefix2+'.txt','w')
        f.write('#centroidInfer.py output:\n')
        f.write('#This catalog is written using inputs from:')
        f.write(irate_file+'\n')
        for i in range(len(arg_list)):
            f.write('#ttype'+str(i)+' = '+arg_list[i]+'\n')

    for i in range(len(subhalo_list)):
        subhalo1= halo1.subhalos[i]

        #here are the intermediate properties that we want to fill our catalog
        #with:
        v3d = subhalo1.props['Velocity']
        subh_center = subhalo1.props['Center']
        subh_center_sigma=subhalo1.props['PositionUncertainty']
        v3d_sigma = subhalo1.props['VelocityUncertainty']
        h0 = halo1.catalog_attrs['h0']
        Ol = halo1.catalog_attrs['Ol']
        Om = halo1.catalog_attrs['Om']
        a = halo1.catalog_attrs['a']
        obs_posx = 0.
        obs_posy = 0.
        obs_posz = 0.
        obs_pos = [obs_posx,obs_posy,obs_posz]
        #have to think about if we want RA increase to +x or -x 
        del_center = subh_center-obs_pos
        del_r = np.linalg.norm(del_center)
        delta = 180./np.pi*np.arcsin(del_center[1]/del_r) 
        #have z=0 to be where RA = 0 
        alpha = 180./np.pi*np.arctan(del_center[0]/del_center[2]) 

        #compute the cosmological redshift from scale factor
        z_cosmo = 1./a - 1. 

        #print 'v3d is {0}'.format(v3d)
        #print 'center is {0}'.format(subh_center)
        #compute radial velocities by first computing the unit vector:
        subh_ucenter = (subh_center - obs_pos)/np.sqrt(np.dot(subh_center-obs_pos,
                                                      subh_center-obs_pos))
        #print 'subh_ucenter is {0}'.format(subh_ucenter)
        # project the 3d velocity unto the radial direction
        v_rad = np.dot(v3d,subh_ucenter)*subh_ucenter
        #divide by the unit vector in order to get the sign right
        v_rad  = (v_rad / subh_ucenter)[0]                                        
        #ignore projection effects from position uncertainty
        v_rad_sigma = subhalo1.props['VelocityUncertainty']
        #----writing out the outputs -----
        if write_out ==True:
            f.write('{0}\t{1}\t{2}\t{3}\n'.format(alpha,delta,v_rad,
                                                     v_rad_sigma))
if write_out==True:
    f.close()
            

#draft of what the script looks like will clean up later



#do error propagation for radial velocity: 

            
