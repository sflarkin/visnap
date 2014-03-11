#! /usr/bin/python
'''Scriptable module with tools to place AHF and Rockstar catalogs into IRATE files'''

import os,sys
from glob import glob
from subprocess import call
from optparse import OptionParser
#import pdb

def cat_and_gen_irate(simdirs, catalog_type, path=None, ahf_catalogs_path = "/analysis/ahf/catalogs/",
                      rock_catalogs_path = "/analysis/Rockstar-0.99.9/output/"):
    '''
    For each dir in simdirs access the catalogs_path, cat the z=0 catalogs from different processes into one catalog,
    and write this catalog into an IRATE file

    Input

     simdirs - a list of simulation directories

     catalog_type - either "ahf" or "rockstar"

     path - if not None simdirs will be searched relative to path, i.e. path/simdirs 

     ahf_catalogs_path - path to ahf catalogs once inside the simdirs

     rock_catalogs_path - path to rockstar catalogs once inside the simdirs
    
    '''

    for thisdir in simdirs:
        simname = thisdir.rpartition('/')[-1]
        snaps = glob(thisdir+"/output/snapshot*")
        if len(snaps) == 0: snaps = glob(thisdir+"/output/"+simname+"*")
        snaps.sort()
        snap = snaps[-1].rpartition('/')[-1]
        snap = snap.rpartition('_')[-1]
        if snap == 'z0': snap = '0'
        if type == 'ahf':
            print '\ncatting AHF_ files for '+simname+' snapshot '+snap
            catdir = path+thisdir+ahf_catalogs_path
            os.system("cat "+catdir+"*.z0.*.AHF_halos > "+catdir+"../"+simname+".z0.AHF_halos")
            os.system("cat "+catdir+"*.z0.*.AHF_profiles > "+catdir+"../"+simname+".z0.AHF_profiles")
            os.system("cat "+catdir+"*.z0.*.AHF_particles > "+catdir+"../"+simname+".z0.AHF_particles")
            print 'generating IRATE file '+simname+'_AHF_irate.hdf5'+' with AHF catalogs'    
            call(["ahf2irate","--profiles","--paramfile="+catdir+simname+".parameter",catdir+"../"+simname+".z0.AHF_",catdir+"../"+simname+"_AHF_irate.hdf5",snap])
        else:
            print '\ncatting rockstar halos_ files for '+simname+' snapshot '+snap
            catdir = path+thisdir+rock_catalogs_path
            os.system("cat "+catdir+"halos_"+snap+".*.ascii > "+catdir+"../"+simname+"_"+snap+".ascii")
            os.system("cat "+catdir+"halos_"+snap+".*.particles > "+catdir+"../"+simname+"_"+snap+".particles")  
            print 'generating IRATE file '+simname+'_rockstar_irate.hdf5'+' with rockstar catalogs'
            call(["rockstar2irate","--param="+catdir+"rockstar.cfg",catdir+"../"+simname+"_"+snap+".ascii",catdir+"../"+simname+"_rockstar_irate.hdf5",snap])
            

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
                      help="Path to IRATE files (default = )")

    (options, args) = parser.parse_args()
    dirs = options.dirs
    path = options.path
    type = options.type

    if dirs == None:
        simdirs = [thisdir for thisdir in os.listdir(path)]
    elif "," in dirs:
        simdirs = dirs.replace(" ","")
        simdirs = simdirs.split(",")
        simdirs.sort()
    else:
        simdirs = [dirs,]

    cat_and_gen_irate(simdirs, type, path)    


    


    
    
    



