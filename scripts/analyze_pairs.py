#! /usr/bin/python
'''Script for analazing halo pairs from the 650 Mpc/h Box'''    

import os, sys, time
from glob import glob
from optparse import OptionParser
import pdb

if __name__ != '__main__':
    print '%s is only script, there is nothing to import from it.' % __name__
    sys.exit()

#Parsing Options

parser = OptionParser()

parser.add_option("-l", "--list", dest="pair_list_file",
                  default = "SIDM_z05_rankedpairs.txt",
                  help="File containig the list of halo pairs")
parser.add_option("-t", "--trees", dest="trees_path",
                  default="catalogs/trees/", 
                  help="Path where trees are located")
parser.add_option("-s", "--snapshot_number", dest="snap",
                  default="194", 
                  help="Snapshot number from which the list of pairs was generated (defaults to 194, which is the z=0.5 snap in the 650 Mpc/h Box run)")
parser.add_option("-L", "--last_snap", dest="Lsnap",
                  default="295", 
                  help="Last snapshot (defaults to 295, which is the last snapshot in the 650 Mpc/h Box run)")

(options, args) = parser.parse_args()
pair_list_file = options.pair_list_file
trees_path = options.trees_path
snap = int(options.snap)
Lsnap = int(options.Lsnap)


#### Start analysis ####
from numpy import *
from visnap.general.halo_track import halo_track, find_tree

#Read list of possible analog pairs 
pair_id, id1, id2, m1, x1, y1, z1, vx1, vy1, vz1, m2, x2, y2, z2, vx2, vy2, vz2, d_3d, weight, phi, theta, dproj, vlos  = genfromtxt(pair_list_file,unpack=True)

#Get the halo1 properties over time
tree_name1, depth_first_ID1 = find_tree(id1[1], snap, Lsnap, trees_path)
halo_past_props1 = halo_track(tree_name1, depth_first_ID1, trees_path)

#Get halo2 properties over time
tree_name2, depth_first_ID2 = find_tree(id1[1], snap, Lsnap, trees_path)
halo_past_props2 = halo_track(tree_name2, depth_first_ID2, trees_path)

#Calculate pair dynamical properties and check if merged

