#! /usr/bin/python
'''Scriptable module with tools to place AHF and Rockstar catalogs into IRATE files'''

import os,sys
from glob import glob
from subprocess import call
from optparse import OptionParser
#import pdb

def cat_and_gen_irate(simdirs, catalog_type, cat_snap_string, snap_num, path=None,
                      ahf_catalogs_path="/analysis/ahf/catalogs/",
                      rock_catalogs_path="/analysis/Rockstar-0.99.9/output/"):
    '''
    For each dir in simdirs access the catalogs_path, cat the catalogs of the given snapshot
    and write this catalog into an IRATE file

    Input

     simdirs - A list of simulation directories

     catalog_type - Either "ahf" or "rockstar"

     cat_snap_string - Snapshot string, the catalogs will be searched as
                       *.cat_snap_string.*.AHF_halos for AHF and halos_cat_snap_string.*.ascii
                       for Rockstar

     snap_num - Snapshot number to use when writing the IRATE files              

     path - If not None simdirs will be searched relative to path, i.e. path/simdirs 

     ahf_catalogs_path - Path to ahf catalogs once inside the simdirs

     rock_catalogs_path - Path to rockstar catalogs once inside the simdirs

    '''
    
    for thisdir in simdirs:
        simname = thisdir.rpartition('/')[-1]
        if catalog_type == 'ahf':
            print '\ncatting AHF_ files for '+simname+' snapshot '+cat_snap_string
            catdir = path+thisdir+'/'+ahf_catalogs_path+'/'
            os.system("cat "+catdir+"*."+cat_snap_string+".*.AHF_halos > "+catdir+"../"+simname+"."+cat_snap_string+".AHF_halos")
            os.system("cat "+catdir+"*."+cat_snap_string+".*.AHF_profiles > "+catdir+"../"+simname+"."+cat_snap_string+".AHF_profiles")
            os.system("cat "+catdir+"*."+cat_snap_string+".*.AHF_particles > "+catdir+"../"+simname+"."+cat_snap_string+".AHF_particles")
            print 'generating IRATE file '+simname+'_AHF_irate.hdf5'+' with AHF catalogs'    
            call(["ahf2irate","--profiles","--paramfile="+catdir+simname+".parameter",catdir+"../"+simname+"."+cat_snap_string+".AHF_",catdir+"../"+simname+"_AHF_irate.hdf5",snap_num])
        else:
            print '\ncatting rockstar halos_ files for '+simname+' snapshot '+cat_snap_string
            catdir = path+thisdir+'/'+rock_catalogs_path+'/'
            os.system("cat "+catdir+"halos_"+cat_snap_string+".*.ascii > "+catdir+"../"+simname+"_"+cat_snap_string+".ascii")
            os.system("cat "+catdir+"halos_"+cat_snap_string+".*.particles > "+catdir+"../"+simname+"_"+cat_snap_string+".particles")  
            print 'generating IRATE file '+simname+'_rockstar_irate.hdf5'+' with rockstar catalogs'
            call(["rockstar2irate","--param="+catdir+"rockstar.cfg",catdir+"../"+simname+"_"+cat_snap_string+".ascii",catdir+"../"+simname+"_rockstar_irate.hdf5",snap_num])
            

if __name__=="__main__":
    
    parser = OptionParser()
    parser.add_option("-d", "--dirs", dest="dirs",
                      default = None,
                      help="Simulation directories to process, e.g. <dir> for a single dir or <dir1,dir2,...> for multiple dirs. If None, will process all directories under 'path' (use -p to set 'path')")
    parser.add_option("-t", "--type", dest="type",
                      default = "rockstar",
                      help="Type of catalogs: ahf or rockstar. (default = rockstar) ")
    parser.add_option("-p", "--path", dest="path",
                      default="", 
                      help="Path to sim directories (default = )")
    parser.add_option("-s", "--snap_string", dest="snap_string",
                      default="",
                      help="Snapshot string, the catalogs will be searched as *.cat_snap_string.*.AHF_halos for AHF and halos_cat_snap_string.*.ascii for Rockstar")
    parser.add_option("-n", "--snap_number", dest="snap_num",
                      default="",
                      help="Snapshot number to use when writing the IRATE files")
    parser.add_option("-a","--ahf_path", dest="ahf_path",
                      default="/analysis/ahf/catalogs/",
                      help="Path to ahf catalogs once inside the simdirs")
    parser.add_option("-r","--rock_path", dest="rock_path",
                      default="/analysis/Rockstar-0.99.9-RC3/output/",
                      help="Path to ahf catalogs once inside the simdirs")

    (options, args) = parser.parse_args()
    dirs = options.dirs
    path = options.path
    type = options.type
    snap_string = options.snap_string
    snap_num = options.snap_num
    ahf_path = options.ahf_path
    rock_path = options.rock_path
    

    if dirs == None:
        simdirs = [thisdir for thisdir in os.listdir(path)]
    elif "," in dirs:
        simdirs = dirs.replace(" ","")
        simdirs = simdirs.split(",")
        simdirs.sort()
    else:
        simdirs = [dirs,]

    cat_and_gen_irate(simdirs, type, snap_string, snap_num, path, ahf_path, rock_path)    


    


    
    
    



