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

from pylab import *
import visnap.general.irate_file_mod as irate_file_mod
import visnap.general.find_halos as find_halos
import visnap.plot.halo_profiles as plot_profiles
import visnap.plot.subhalos as plot_subhalos

# Initialize halo objects

halos = []   # List of halo objects
for irate_file in irate_files:
    snap_names, redshifts, scales  = irate_file_mod.snapshot_times(irate_file)
    z0snap = snap_names[scales == 1][0]
    mostpart_halo = find_halos.find_mostpart_halo(irate_file, z0snap, 'AHF')
    halos.append(mostpart_halo)


'''
# Set legend names
legendNames = [r'$\mathrm{CDM}$', r'$\mathrm{SIDM} \ \sigma/m = 0.1 \ \mathrm{cm^2/g}$',
               r'$\mathrm{SIDM} \ \sigma/m = 1 \ \mathrm{cm^2/g}$',
               r'$\mathrm{WDM} \ m = 15 \ \mathrm{kev}$', r'$\mathrm{WDM} \ m = 6 \ \mathrm{kev}$',
               r'$\mathrm{DDM} \ \tau = 10  \ \mathrm{gyr}, \ v_k = 45  \ \mathrm{km/s}$', 
               r'$\mathrm{DDM} \ \tau = 1  \ \mathrm{gyr}, \ v_k = 20  \ \mathrm{km/s}$']
'''

'''
# Density profiles
fig, ax, lines, legends = plot_profiles.density_profiles(halos)
ax.set_xlim(1,200)
ax.set_ylim(5e-6,.2)
#legend = ax.legend(lines,legendNames,loc=3,prop=dict(size=0.8*25,),labelspacing = 0.1)
#legend.draw_frame(False)
ax.set_yticks(logspace(-5,-1,5))
ax.set_yticklabels([r'$-5$',r'$-4$',r'$-3$',r'$-2$',r'$-1$'],size=25)
ax.set_xticks([1,5,10,50,100,200])
ax.set_xticklabels([r'$1$',r'$5$',r'$10$',r'$50$',r'$100$',r'$200$'],size=25)
show()

# Subhalo function
fig, ax, lines, legends = plot_subhalos.subhalo_functions(halos, rcut_factor = 1, vcut_factor = 0.05)
ax.set_xlim(8,40)
ax.set_ylim(0,150)
#legend = ax.legend(lines,legendNames,loc=1,prop=dict(size=0.7*25,),labelspacing = 0.1)
#legend.draw_frame(False)
ax.set_yticks([1,5,10,50,100])
ax.set_yticklabels([r'$1$',r'$5$',r'$10$',r'$50$',r'$100$'],size=25)
ax.set_xticks([8,10,20,30,40])
ax.set_xticklabels([r'$8$',r'$10$',r'$20$',r'$30$',r'$40$'],size=25)
show()
'''

# Subhalo 2Dfunction
fig, ax, lines, legends = plot_subhalos.subhalo_2DRfunctions(halos, projection_axis='minor',
                                                             rcut_factor = 2,
                                                             vcut_factor = 0.05,
                                                             showme=0)
ax.set_xlim(0.01,2)
ax.set_ylim(0.8,200)
#legend = ax.legend(lines,legendNames,loc=1,prop=dict(size=0.7*25,),labelspacing = 0.1)
#legend.draw_frame(False)
ax.set_yticks([1,5,10,50,100])
ax.set_yticklabels([r'$1$',r'$5$',r'$10$',r'$50$',r'$100$'],size=25)
ax.set_xticks([0.01,0.05,0.1,0.5,1,2])
ax.set_xticklabels([r'$0.01$',r'$0.05$',r'$0.1$',r'$0.5$',r'$1$',r'$2$'],size=25)
show()









