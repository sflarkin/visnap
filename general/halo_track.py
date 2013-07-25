'''Module containing tools for tracking halos over time'''

import os, subprocess, sys, glob
import h5py
from numpy import *
from multiprocessing import Process, Queue, cpu_count
from visnap.general.rockstar_trees2irate import trees2irate
import pdb

def find_tree(orig_halo_id, orig_snapshot, Lsnap, trees_path,
              trees_file_base='tree', ncpus='all'):
    '''
    Find the tree containing a given halo. For now this routine only works for
    rockstar halo catalogs and trees

    Input:

     orig_halo_id - Rockstar halo ID in catalog of orig_snapshot 

     orig_snapshot - Snapshot from which orig_halo_id was originated

     Lsnap - Last snapshot number for this run. Used to figure out how
             many snapshots the rockstar tree generator skipped due to 
             not enough halos.

     trees_path - The path to the merger trees files

     trees_file_base - The base name of the  merger trees files

     ncpus - If 'all' all the avilable cores in the node will be used, else set
             to the number of cores that you want to use

    Output:

     found_tree - The h5py key name of the found tree. Will be empty if a tree
                  wasn't found

     depth_first_id - Depth first ID of the halo found in the tree 

    '''
    
    #First we want to see if the merger trees IRATE file exists, if not we generate it.
    if not os.path.exists(trees_path+trees_file_base+'.hdf5'):
        print "No IRATE file %s containing trees found, now generating" % (trees_path+trees_file_base+'.hdf5')         
        Ntree_files = len(glob.glob(trees_path+trees_file_base+'_*.dat'))
        if Ntree_files > 1:
            subprocess.call(['cat '+trees_path+trees_file_base+'_*.dat'+' > '+trees_path+trees_file_base+'.dat'], shell=True)
            trees2irate(trees_path+trees_file_base+'.dat',trees_path+trees_file_base+'.hdf5')
            subprocess.call(['rm '+trees_path+trees_file_base+'.dat'], shell=True)
        elif Ntree_files:
            trees2irate(trees_path+trees_file_base+'.dat',trees_path+trees_file_base+'.hdf5')
        else:
            print "The Rockstar trees file %s was not found" % (trees_path+trees_file_base+'_*.dat')
            sys.exit()

    #Now we want to look for the tree containing this halo. We do this in parallel using ncores.
    
    if ncpus == 'all':
        ncores = cpu_count()
    else:
        ncores = int(ncpus)
    
    print 'Starting tree search with %d cpus\n' % ncores

    def find_halo_intree(orig_halo_id,orig_snapshot,Lsnap,in_que,out_que):
        '''
        Search for halo in the trees

        Kernel function to be run in multiple processes
        '''
        f = h5py.File(trees_path+trees_file_base+'.hdf5')
        trees = f['RockstarMergerTrees']
        for keys in iter(in_que.get,'STOP'):
            for key in keys:
                tree = trees[key]
                offset = Lsnap-tree['Snap_num'][...].max()
                snap = orig_snapshot-offset
                arg = argwhere((tree['Orig_halo_ID'][...] == orig_halo_id) & (tree['Snap_num'][...] == snap))
                if len(arg): 
                    print 'Halo %d found in %s' % (orig_halo_id,key)
                    out_que.put((key, tree['Depth_first_ID'][...][arg[0,0]]))
                    break 
           
    # start ncores process calling the above kernel function 
    in_que = Queue()
    out_que = Queue()
    jobs = []
    for i in range(ncores):
        job = Process(target = find_halo_intree, args = (orig_halo_id,orig_snapshot,Lsnap,in_que,out_que))
        job.start()
        jobs.append(job)

    # fill up the in_que with chunks of trees, we actually pass only the tree names.  
    f = h5py.File(trees_path+trees_file_base+'.hdf5')
    trees = f['RockstarMergerTrees']
    Ntrees = len(trees)    
    trees_chunk = []
    chunk_size = Ntrees/(ncores*10)
    for i,tree in enumerate(trees.iterkeys()):
        if out_que.qsize() > 0:
           #the tree was already found 
           break
        else:
            trees_chunk.append(tree)
            if (len(trees_chunk) == chunk_size) or (i == Ntrees-1):
                in_que.put(trees_chunk)
                trees_chunk = []                         
      
    # we are done send the stop signal to all processes, and wait for them to terminate and join.
    for i in range(ncores): in_que.put('STOP')
    for job in jobs: job.join()
    
    # If the tree was found the out_que should have only one tree, but we check for
    # for multiple trees returned just in case    
    found_tree = [out_que.get() for i in range(out_que.qsize())]
    if len(found_tree) == 0:
        print 'A tree containing halo %d in snaphshot %d was not found'%(orig_halo_id,orig_snapshot)                 
    elif len(found_tree) > 1:
        print 'WARNING: Multiple trees found to contain halo %d in snapshot %d'%(orig_halo_id,orig_snapshot) 
        print 'These are the trees found: ',[ft[0] for ft in found_tree] 
        print 'I will return the first one'
        found_tree = found_tree[0]
    else:
        found_tree = found_tree[0]
        print 'Succesfully found halo %d in tree %s from snapshot %d'%(orig_halo_id,found_tree[0],orig_snapshot)
    #close queues
    in_que.close()
    out_que.close()
    
    #Now found_tree[0] is the tree we want and found_tree[1] is the
    #depth_first_ID of the halo found in the tree 
    return found_tree
        

def halo_track(tree_name, halo_depth_first_ID, trees_path,
               trees_file_base='tree', halo_track='MMP'):

    '''
    Trace halo properties back in time. For now this routine only works for
    rockstar halo catalogs and trees

    Input:
 
     tree_name - The h5py key name of the tree containig the halo    
     
     halo_depth_first_id - The rockstar Depth_first_ID of the halo we want to track
     
     trees_path - The path to the merger trees files
     
     trees_file_base - the base of the  merger trees files
     
     halo_track - 'MMP' to follow the most massive progenitor track, or 'MMA' to follow the 
                  track with most mass above ('MMA not implemented yet').

    Output:
     
     halo_past_props - A dictionary containing the halo properties over time,
                       all halo properties in the rockstar tree are included 
                       with the same key names. 
    '''
   
    f = h5py.File(trees_path+trees_file_base+'.hdf5')
    trees = f['RockstarMergerTrees']
    tree = trees[tree_name]

    df_id = tree['Depth_first_ID'][...]
    scale = tree['scale'][...]
    halo_past_props = {}

    # We first start the dictionary with the properties at the current time of the halo
    for key in tree.keys():
        if key == 'T': continue
        halo_past_props[key] = tree[key][...][df_id == halo_depth_first_ID]
    
    # Then we use the depth_first_ID to go back in time followint the MMP and record the halo properties  
    # at each time step    
    old_scale = halo_past_props['scale'][-1]
    new_df_id = halo_depth_first_ID + 1    
    new_scale = scale[df_id == new_df_id][0]   
    max_df_id = df_id.max()
    while((old_scale > new_scale) and (new_df_id <= max_df_id)):
        for key in tree.keys():
            if key == 'T': continue
            halo_past_props[key] = append(halo_past_props[key], tree[key][...][df_id == new_df_id])
        old_scale = new_scale
        new_df_id = new_df_id + 1
        new_scale = scale[df_id == new_df_id ]  

    # And we got it, now halo_past_props contains the properties of the halos that form the MMP branch of our halo
    return halo_past_props



