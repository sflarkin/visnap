#! /usr/bin/python
'''Template useful to start new analysis scripts utilizing VISNAP modules'''    

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

#Do whatever you want below

#example
import visnap.general.translate_filename as translate_filename
import visnap.general.find_halos as find_halos
import visnap.plot.halo_profiles as plot_profiles

halos = []
for irate_file in irate_files:
    simprops = translate_filename.hyades(irate_file)
    zoom_id = simprops['zoom_id']
    zoom_halo, dist = find_halos.find_zoom_halo('../650Box_clusters_zooms.txt', zoom_id,
                                                irate_file, 'Snapshot00148', 'Rockstar')
    halos.append(zoom_halo)
    
#fig, ax, lines, legends = plot_profiles.density_profiles(halos)
 


