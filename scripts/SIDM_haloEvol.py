#! /usr/bin/python
'''
Script to analyse halo density profiles and rotation curves over time, 
and to compare them between SIDM and CDM halos.
'''    

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
                  "file and <path1/file1,path2/file2,...> or <'path1/file1, "\
                  "path2/file2, ...'>  for multiple files. Leave blank to read "\
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
elif "," in irate_files:
    irate_files = irate_files.split(",")
    irate_files = [irate_file.strip(' ') for irate_file in irate_files]
    irate_files.sort()    
else:
    irate_files = [irate_files,]

print "\nThe following %d IRATE files were parsed for analysis:\n"%len(irate_files),irate_files[:]
print '' 

#Do whatever you want below

#Import stuff we will need
import h5py
from numpy import *
import visnap.plot
import visnap.general.find_halos as find_halos
import visnap.plot.halo_profiles as plot_profiles

#Get halos for analysis
halos = []
for irate_file in irate_files:
    for snap in range(10, 16):
        halo = find_halos.find_mostpart_halo(irate_file, snap, 'Rockstar')
        halos.append(halo)

#Now halos is a list of halo objects containing the halo with the most number 
#of particles for snapshots 1-15 in all the given IRATE files

#Assign a label, color and marker to each of these halos (for plotting)  
labels = []
markers = []
for halo in halos:
    dmName = halo.sim_props['dmName']
    if dmName == 'CDM': marker = '-'
    else: marker = '--'
    markers.append(marker)

    redshift = 1/halo.catalog_attrs['a'] - 1
    Vmax = halo.props['Vmax']
    Mvir = halo.props['Mvir']    
    labels.append(dmName+'\_'+str(round(redshift,2))+'\_'+str(round(Vmax))[:-2])
    
colors_list = visnap.plot.colors_list
nhalos = len(halos)
colors = concatenate((colors_list[0:nhalos/2], colors_list[0:nhalos/2]))

#Select which halos to plot

fig, ax, lines, legends = plot_profiles.density_profiles(halos, colors=colors, legendnames=labels)
fig, ax, lines, legends = plot_profiles.circular_velocities(halos, colors=colors,legendnames=labels)


